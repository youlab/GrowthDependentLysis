import numpyro
import numpyro.distributions as dist
from numpyro.infer import NUTS, MCMC
import numpy as np
from jax import numpy as jnp
import jax.random as random
import pickle


class BayesLinearRegress:
    def __init__(self, X, y,
        use_intercept=True, rescale_data=True,
        MCMC_samples=2000, MCMC_warmup=2000,
        param_hyperprior = 1e0,
        measurement_noise_hyperprior = 1e0,
    ):
        """
        Args:
            X: Data (n_samples, n_features)
            y: Target (n_samples,)
            use_intercept: Whether to add an intercept term to the design matrix
            rescale_data: Whether to normalize the data
            MCMC_samples: Number of posterior samples to draw
            MCMC_warmup: Number of warmup (burn-in) steps for NUTS
            param_hyperprior: Scale of HalfCauchy prior on beta's std (tau)
            measurement_noise_hyperprior: Scale of HalfCauchy prior on noise std (sigma)
        """
        self.use_intercept = use_intercept
        self.rescale_data = rescale_data
        self.X = X # Data
        self.y = y # Target

        # Add intercept if needed
        if use_intercept:
            self.X = np.column_stack((np.ones(len(X)), X))

        # Normalize data where X = (X - X_mean) / X_std
        self.X_mean = np.mean(self.X, axis=0)
        self.X_std = np.std(self.X, axis=0) + 1e-6
        self.y_mean = np.mean(self.y)
        self.y_std = np.std(self.y) + 1e-6
        if rescale_data:
            self.X_rescaled = (self.X - self.X_mean) / self.X_std
            self.y_rescaled = (self.y - self.y_mean) / self.y_std
        else:
            self.X_rescaled = self.X
            self.y_rescaled = self.y

        # Set the hyperparameters and model parameters
        self.MCMC_samples = MCMC_samples
        self.MCMC_warmup = MCMC_warmup
        self.param_hyperprior = param_hyperprior
        self.measurement_noise_hyperprior = measurement_noise_hyperprior
        self.beta = None  # MCMC samples, populated after fit()
        self.sigma = None # MCMC samples, populated after fit()
        self.tau = None   # MCMC samples, populated after fit()
        self.num_features = self.X.shape[1]  # includes intercept column if added

    def fit(self, rng_seed=0):
        """
        Run NUTS MCMC to sample from the posterior distribution.

        NumPyro workflow:
        1. NUTS builds a Hamiltonian Monte Carlo kernel that uses JAX autodiff
           to compute gradients of the log-posterior, and automatically tunes
           step size and trajectory length.
        2. MCMC runs warmup steps (adapting sampler params) then draws
           posterior samples.
        3. Samples for beta, sigma, tau are extracted and stored as numpy arrays.

        Args:
            rng_seed: Random seed for JAX PRNG (for reproducibility)
        """
        kernel = NUTS(self.MCMC_model)
        mcmc = MCMC(kernel, num_warmup=self.MCMC_warmup, num_samples=self.MCMC_samples)
        mcmc.run(
            random.PRNGKey(rng_seed),
            X=jnp.array(self.X_rescaled),
            y=jnp.array(self.y_rescaled),
        )
        samples = mcmc.get_samples()
        self.beta = np.array(samples["beta"])    # (MCMC_samples, num_features)
        self.sigma = np.array(samples["sigma"])  # (MCMC_samples,)
        self.tau = np.array(samples["tau"])       # (MCMC_samples,)

    def create_kernel(self):
        """Return a NUTS kernel bound to this instance's MCMC model.

        The kernel can be reused across multiple MCMC runs (JIT-compiled once).
        """
        return NUTS(self.MCMC_model)

    def run_single(self, kernel, X_raw, y_raw, rng_key, num_warmup=None, num_samples=1):
        """Run MCMC with a pre-built kernel on new data, returning raw samples.

        Builds the design matrix, rescales using stored stats, and runs MCMC.
        Useful for composition sampling where the same model structure is
        re-fitted many times on different data realizations.

        Args:
            kernel: NUTS kernel from create_kernel()
            X_raw: (n, n_features) raw feature data (without intercept)
            y_raw: (n,) raw target data
            rng_key: JAX PRNGKey
            num_warmup: warmup steps (defaults to self.MCMC_warmup)
            num_samples: number of posterior samples to draw (default 1)

        Returns:
            dict of raw MCMC samples in rescaled space (beta, sigma, tau)
        """
        X = np.array(X_raw)
        y = np.array(y_raw)

        if self.use_intercept:
            X = np.column_stack((np.ones(len(X)), X))

        if self.rescale_data:
            X = (X - self.X_mean) / self.X_std
            y = (y - self.y_mean) / self.y_std

        mcmc = MCMC(
            kernel,
            num_warmup=num_warmup or self.MCMC_warmup,
            num_samples=num_samples,
            progress_bar=False,
        )
        mcmc.run(rng_key, X=jnp.array(X), y=jnp.array(y))
        return mcmc.get_samples()

    def predict(self, X):
        """
        Draw from the posterior predictive distribution for new data.

        For each MCMC sample (beta_i, sigma_i), draws:
            y_i ~ N(X @ beta_i, sigma_i)
        then inverse-transforms back to the original y space.

        Args:
            X: New data (n_samples, n_features) in original (unscaled) space
        Returns:
            y_pred: Posterior predictive samples (n_samples, MCMC_samples)
                    in original space.
        """
        if self.beta is None:
            raise RuntimeError("Model has not been fitted. Call fit() first.")

        if self.use_intercept:
            X = np.column_stack((np.ones(len(X)), X))

        if self.rescale_data:
            X = (X - self.X_mean) / self.X_std

        # (n_samples, num_features) @ (num_features, MCMC_samples) -> (n_samples, MCMC_samples)
        mu_rescaled = X @ self.beta.T

        # Sample from the likelihood: y ~ N(mu, sigma) for each posterior draw
        # sigma shape: (MCMC_samples,) -> broadcast with (n_samples, MCMC_samples)
        noise = np.random.randn(*mu_rescaled.shape) * self.sigma[np.newaxis, :]
        y_rescaled = mu_rescaled + noise

        if self.rescale_data:
            y_pred = y_rescaled * self.y_std + self.y_mean
        else:
            y_pred = y_rescaled

        return y_pred

    def MCMC_model(self, X=None, y=None):
        """
        NumPyro probabilistic model for Bayesian linear regression.

        This function defines the generative model that NumPyro uses for MCMC.
        Each numpyro.sample() call declares a random variable in the model.
        The obs=y argument on the last sample tells NumPyro that y is observed
        data, so MCMC infers the posterior of (tau, sigma, beta) conditioned
        on the observed y.

        Generative process:
            tau ~ HalfCauchy(param_hyperprior)
            sigma ~ HalfCauchy(measurement_noise_hyperprior)
            beta ~ N(0, tau^2 * I)
            y | X, beta, sigma ~ N(X @ beta, sigma)

        Args:
            X: Design matrix, passed automatically by MCMC.run()
            y: Observed targets, passed automatically by MCMC.run()
        """
        # Sample the parameter hyperprior
        tau = numpyro.sample("tau", dist.HalfCauchy(scale=self.param_hyperprior))
        sigma = numpyro.sample("sigma", dist.HalfCauchy(scale=self.measurement_noise_hyperprior))

        # Sample the parameters
        beta = numpyro.sample(
            "beta", dist.MultivariateNormal(
                loc=jnp.zeros(self.num_features),
                covariance_matrix=tau**2 * jnp.eye(self.num_features)
            )
        )

        # obs=y conditions the model on observed data:
        # instead of sampling y from the prior, NumPyro evaluates the
        # log-likelihood of the observed y given X @ beta and sigma.
        numpyro.sample("y", dist.Normal(X @ beta, scale=sigma), obs=y)

    def save(self, filepath):
        """
        Save model config, normalization stats, and MCMC samples to file.
        Only stores what is needed for predict() — not the training data.

        Args:
            filepath: Path to save the model (e.g. "model.pkl")
        """
        state = {
            "use_intercept": self.use_intercept,
            "rescale_data": self.rescale_data,
            "X_mean": self.X_mean,
            "X_std": self.X_std,
            "y_mean": self.y_mean,
            "y_std": self.y_std,
            "MCMC_samples": self.MCMC_samples,
            "MCMC_warmup": self.MCMC_warmup,
            "param_hyperprior": self.param_hyperprior,
            "measurement_noise_hyperprior": self.measurement_noise_hyperprior,
            "num_features": self.num_features,
            "beta": self.beta,
            "sigma": self.sigma,
            "tau": self.tau,
        }
        with open(filepath, "wb") as f:
            pickle.dump(state, f)

    @classmethod
    def load(cls, filepath):
        """
        Load a saved model from file. The returned instance can call predict()
        without needing the original training data.

        Args:
            filepath: Path to the saved model file
        Returns:
            BayesLinearRegress instance with loaded state
        """
        with open(filepath, "rb") as f:
            state = pickle.load(f)

        model = cls.__new__(cls)
        model.use_intercept = state["use_intercept"]
        model.rescale_data = state["rescale_data"]
        model.X_mean = state["X_mean"]
        model.X_std = state["X_std"]
        model.y_mean = state["y_mean"]
        model.y_std = state["y_std"]
        model.MCMC_samples = state["MCMC_samples"]
        model.MCMC_warmup = state["MCMC_warmup"]
        model.param_hyperprior = state["param_hyperprior"]
        model.measurement_noise_hyperprior = state["measurement_noise_hyperprior"]
        model.num_features = state["num_features"]
        model.beta = state["beta"]
        model.sigma = state["sigma"]
        model.tau = state["tau"]
        model.X = None
        model.y = None
        model.X_rescaled = None
        model.y_rescaled = None
        return model


if __name__ == "__main__":
    import os
    import matplotlib.pyplot as plt

    # --- Synthetic 1D data: y = 2*x + 5 + noise ---
    np.random.seed(42)
    n = 50
    X_train = np.random.uniform(-3, 3, (n, 1))
    true_slope = 2.0
    true_intercept = 5.0
    noise_std = 2.8
    y_train = true_slope * X_train[:, 0] + true_intercept + noise_std * np.random.randn(n)

    # Fit
    model = BayesLinearRegress(X_train, y_train, MCMC_samples=500, MCMC_warmup=500)
    model.fit(rng_seed=0)

    # Dense grid for plotting
    x_grid = np.linspace(-4, 4, 200).reshape(-1, 1)
    y_pred = model.predict(x_grid)  # (200, 500)
    y_mean = y_pred.mean(axis=1)
    y_lo = np.percentile(y_pred, 2.5, axis=1)
    y_hi = np.percentile(y_pred, 97.5, axis=1)
    y_true = true_slope * x_grid[:, 0] + true_intercept

    print("=== Bayesian Linear Regression Test ===")
    print(f"True: slope={true_slope}, intercept={true_intercept}, noise_std={noise_std}")
    print(f"Posterior sigma mean: {model.sigma.mean():.4f}")

    # --- Plot ---
    fig, axes = plt.subplots(1, 4, figsize=(18, 4))

    # 1) Regression fit with 95% credible band
    ax = axes[0]
    ax.scatter(X_train[:, 0], y_train, s=20, alpha=0.6, color="black", label="Training data")
    ax.plot(x_grid[:, 0], y_true, "--", color="red", linewidth=1.5, label="True line")
    ax.plot(x_grid[:, 0], y_mean, color="blue", linewidth=1.5, label="Posterior mean")
    ax.fill_between(x_grid[:, 0], y_lo, y_hi, alpha=0.25, color="blue", label="95% credible interval")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Bayesian Linear Regression Fit")
    ax.legend(fontsize=8)

    # 2) Posterior draws — plot a handful of sampled lines
    ax = axes[1]
    ax.scatter(X_train[:, 0], y_train, s=20, alpha=0.4, color="black")
    for j in range(min(50, y_pred.shape[1])):
        ax.plot(x_grid[:, 0], y_pred[:, j], color="blue", alpha=0.06)
    ax.plot(x_grid[:, 0], y_true, "--", color="red", linewidth=1.5, label="True line")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Posterior Draws")
    ax.legend(fontsize=8)

    # 3) Posterior histogram for sigma
    ax = axes[2]
    ax.hist(model.sigma, bins=30, alpha=0.6, color="steelblue")
    ax.axvline(noise_std, color="red", linestyle="--", linewidth=1.5, label=f"True = {noise_std}")
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel("Count")
    ax.set_title(r"Posterior $\sigma$ (noise std)")
    ax.legend(fontsize=8)

    # 4) Posterior histogram for tau
    ax = axes[3]
    ax.hist(model.tau, bins=30, alpha=0.6, color="coral")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("Count")
    ax.set_title(r"Posterior $\tau$ (prior std)")

    plt.tight_layout()
    plt.savefig("bayesian_linreg_test.pdf")
    plt.show()
    print("Plot saved to bayesian_linreg_test.pdf")

    # Test save/load roundtrip (compare stored samples, not predictions,
    # since predict() draws fresh noise each call)
    save_path = "test_model.pkl"
    model.save(save_path)
    loaded_model = BayesLinearRegress.load(save_path)
    params_match = (
        np.allclose(model.beta, loaded_model.beta)
        and np.allclose(model.sigma, loaded_model.sigma)
        and np.allclose(model.tau, loaded_model.tau)
    )
    print(f"Save/load roundtrip: parameters match = {params_match}")
    os.remove(save_path)
