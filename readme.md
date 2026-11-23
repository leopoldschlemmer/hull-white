# Hull–White model fitted to Swiss government bonds

The calculations below derive the short-rate simulation and the function $\phi$ used in [calibration.ipynb](calibration.ipynb). The notebook constructs a discount curve from SNB Swiss Confederation spot yields and chooses $\phi$ to match that initial curve.

## Run the notebook

Use Python 3.11 or newer. From the repository directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m jupyterlab calibration.ipynb
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. Select the environment's Python kernel, then run all notebook cells in order.

The repository includes the raw [SNB bond-yield dataset](https://data.snb.ch/en/topics/ziredev/cube/rendoblid) and its [source metadata](data/source_metadata.json), so the calculations run offline after installing dependencies. The default example uses the curve from 4 January 1988, 5,000 paths, and a fixed random seed. Change `valuation_date`, `a`, and `sigma` in the simulation cell to explore another case.

## 1. Short-rate model

Under the risk-neutral measure $\mathbb{Q}$, write

$$
r_t = \phi(t) + X_t,
\qquad
dX_t = -aX_t\,dt + \sigma\,dW_t^{\mathbb{Q}},
\qquad X_0 = 0.
$$

Here $\phi$ is deterministic, $a>0$ is the mean-reversion parameter, and $\sigma\geq 0$ is the volatility. We first solve for $X$, then use bond prices to determine $\phi$.

## 2. Distribution of X

Multiplying by the integrating factor $e^{at}$ gives

$$
\begin{aligned}
d(e^{at}X_t)
&= ae^{at}X_t\,dt + e^{at}\,dX_t \\
&= ae^{at}X_t\,dt - ae^{at}X_t\,dt
   + \sigma e^{at}\,dW_t^{\mathbb{Q}} \\
&= \sigma e^{at}\,dW_t^{\mathbb{Q}}.
\end{aligned}
$$

Since $X_0=0$,

$$
e^{at}X_t = \sigma\int_0^t e^{as}\,dW_s^{\mathbb{Q}},
\qquad
X_t = \sigma\int_0^t e^{-a(t-s)}\,dW_s^{\mathbb{Q}}.
$$

The integrand is deterministic, so $X_t$ is Gaussian with mean zero. By Itô's isometry,

$$
\begin{aligned}
\operatorname{Var}^{\mathbb{Q}}(X_t)
&= \sigma^2\int_0^t e^{-2a(t-s)}\,ds \\
&= \sigma^2 e^{-2at}\frac{e^{2at}-1}{2a} \\
&= \frac{\sigma^2}{2a}(1-e^{-2at}).
\end{aligned}
$$

Therefore,

$$
X_t\sim\mathcal{N}\left(0,\frac{\sigma^2}{2a}(1-e^{-2at})\right).
$$

## 3. Conditional distribution and simulation

Integrating between $t$ and $t+h$ instead gives

$$
e^{a(t+h)}X_{t+h}-e^{at}X_t
=\sigma\int_t^{t+h}e^{as}\,dW_s^{\mathbb{Q}},
$$

and hence

$$
X_{t+h}
=e^{-ah}X_t
+\sigma\int_t^{t+h}e^{-a(t+h-s)}\,dW_s^{\mathbb{Q}}.
$$

The second term is independent of the information at time $t$. Thus,

$$
X_{t+h}\mid X_t=x
\sim\mathcal{N}\left(e^{-ah}x,\frac{\sigma^2}{2a}(1-e^{-2ah})\right).
$$

Starting at $X_0=0$, we simulate on a grid using independent standard normal variables $Z_i$:

$$
X_{t_{i+1}}
=e^{-ah_i}X_{t_i}
+\sqrt{\frac{\sigma^2}{2a}(1-e^{-2ah_i})}\,Z_i,
\qquad h_i=t_{i+1}-t_i.
$$

These are exact transitions for $X$. The short-rate paths follow from $r_{t_i}=\phi(t_i)+X_{t_i}$.

## 4. Bond prices

The price of a zero-coupon bond paying one unit at maturity $T$ is

$$
B(0,T)=\mathbb{E}^{\mathbb{Q}}\left[
\exp\left(-\int_0^T r_s\,ds\right)
\right].
$$

Using $X_t=r_t-\phi(t)$, the equation for $X$ becomes

$$
dX_t=-ar_t\,dt+a\phi(t)\,dt+\sigma\,dW_t^{\mathbb{Q}}.
$$

Integrating and rearranging, with $X_0=W_0^{\mathbb{Q}}=0$,

$$
\begin{aligned}
\int_0^T r_s\,ds
&=\int_0^T\phi(s)\,ds
  +\frac{\sigma}{a}W_T^{\mathbb{Q}}-\frac{X_T}{a} \\
&=\int_0^T\phi(s)\,ds
  +\frac{\sigma}{a}\int_0^T
     (1-e^{-a(T-s)})\,dW_s^{\mathbb{Q}}.
\end{aligned}
$$

This is Gaussian with mean

$$
\mu_T=\int_0^T\phi(s)\,ds
$$

and variance

$$
v_T=\frac{\sigma^2}{a^2}\int_0^T(1-e^{-a(T-s)})^2\,ds.
$$

For a Gaussian variable $Y\sim\mathcal{N}(\mu,v)$, we have $\mathbb{E}[e^{-Y}]=e^{-\mu+v/2}$. Therefore,

$$
B(0,T)=\exp\left(
\frac{\sigma^2}{2a^2}\int_0^T(1-e^{-a(T-s)})^2\,ds
-\int_0^T\phi(s)\,ds
\right).
$$

## 5. Choosing phi to fit the initial curve

Taking logarithms and rearranging gives

$$
\int_0^T\phi(s)\,ds
=\frac{\sigma^2}{2a^2}\int_0^T(1-e^{-a(T-s)})^2\,ds
-\log B(0,T).
$$

The substitution $u=T-s$ shows that

$$
\int_0^T(1-e^{-a(T-s)})^2\,ds
=\int_0^T(1-e^{-au})^2\,du.
$$

Differentiating with respect to maturity $T$, we obtain

$$
\boxed{
\phi(T)
=-\frac{\partial_T B(0,T)}{B(0,T)}
+\frac{\sigma^2}{2a^2}(1-e^{-aT})^2
}.
$$



## 6. Connection to the notebook

The notebook converts spot yields into discount factors and linearly interpolates the corresponding bond prices. It uses the interval slope for $\partial_T B(0,T)$ between observations. At interior grid points, where the interpolated curve has a kink, it uses the average of the two adjacent slopes.

All bond prices are for a unit payoff, with $B(0,0)=1$. The interval from zero to the first observed maturity is also interpolated linearly.

After simulating the short rates with the same $a$ and $\sigma$ used in $\phi$, bond prices are checked by Monte Carlo:

$$
\widehat{B}(0,T)
=\frac{1}{N}\sum_{j=1}^N
\exp\left(-\int_0^T r_s^{(j)}\,ds\right).
$$

The time integrals are approximated using the trapezoidal rule. Although the transitions for $X$ are exact, the numerical bond prices still have time-integration error and Monte Carlo sampling error.


## Acknowledgements
AI assistance: The mathematical derivations are based on my own notes. I used ChatGPT to help turn those notes into this write-up and to assist with the Python implementation.