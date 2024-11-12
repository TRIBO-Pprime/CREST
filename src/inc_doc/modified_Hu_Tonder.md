
##Principle

The well-known Hu and Tonder **[8]** procedure is reminded hereafter, with a slight improvement proposed by Bakolas **[3]**.

###Step 1
*Determination of the desired spatial characteristics (acf), and the four statistical moments \(\mu_i\) (mean \(\mu\), variance \(\sigma^2\), skewness Sk and kurtosis Ku)*

Unlike *Skz* and *Kuz* which are chosen by the user in order to be representative of the desired worn surface, the prescribed mean \(\mu\) is set to 0 and the standard deviation \(\sigma\) to 1.
However, at the end of the surface generation process, these two moments can be changed with no effect on the third and fourth moments, by scaling and shifting the final surface heights.

###Step 2
*Determination of the digital filter H.*

The use of a digital filter writes \(z=h \otimes \eta\) where \(z\) is the final surface height, \(\eta\) a random white noise and \(h\) the digital filter.
In the frequency space it becomes \(Z=H \cdot A\), then, \(Z\bar{Z}=H\bar{H}A\bar{A}\) which is simply written as \(Z^2=H^2 A^2\). As \(Z^2=\text{FFT}(acf)\),
and because \(h\) is symmetric regarding \(x\) and \(y\), \(h=\text{FFT}^{-1} \left[ \sqrt{\left|\text{FFT}(acf_z)\right|} \right]/|A|\).
If \(\eta\) is white noise, \(|A|\) becomes a constant and can be ignored.

###Step 3
*Determination of the starting values \(Sk_\eta\) and \(Ku_\eta\), derived from the desired values \(Sk_z\) and \(Ku_z\)*

A starting random set \(\eta\) is generated. \(\eta\) is considered as nearly white noise, centered and scaled: \(\mu_\eta=0\), \(\sigma_\eta=1\).
If a digital filter \(h\) is applied to \(\eta\), the \(z\) resulting statistical moments are modified except for the mean.

\begin{align*}
	\mu_z      &   =    0                                                                                    \\
	\sigma_z^2 &\approx \sum_{k=1}^{n} h_k^2                                                                 \\
	Sk_z       &\approx \dfrac{\sum_{k=1}^{n} h_k^3}{\left[\sum_{k=1}^{n} h_k^2 \right]^\frac{3}{2}} Sk_\eta \\
	Ku_z -3    &\approx \dfrac{\sum_{k=1}^{n} h_k^4}{\left[\sum_{k=1}^{n} h_k^2 \right]^2          } (Ku_\eta -3)
\end{align*}

In this work, the authors propose the use of an alternative set of equations rather than the classical one above.
It is simpler and uses the zero-mean property of  \(\eta\). Indeed  \(h\) can be written \(h=\bar{h} + \mu_h\) where  \(\bar{h}\) is zero-mean,
then \(z=(\bar{h} + \mu_h) \otimes \eta = \bar{h} \otimes \eta\):

\begin{align*}
	\mu_z      &  =  0                                            \\
	\sigma_z^2 &\approx n \sigma_h^2                              \\
	Sk_z       &\approx \dfrac{1}{\sqrt{n}} Sk_h Sk_\eta          \\
	Ku_z -3    &\approx \frac{1}{n}        (Ku_h -3)(Ku_\eta -3)
\end{align*}

###Step 4
*Johnson's transformation of \(\eta\) to match \(Sk_\eta\) and \(Ku_\eta\)*

The Johnson’s translation system is usually utilized to transform non-Gaussian data sets into Gaussian data sets.
In the present case, the system is reversed in order to obtain data sets \(\eta\) with prescribed moments, starting from Gaussian noise \(\phi\).
The translator system is composed of three functions reminded hereafter:

\begin{align*}
	S_U : \eta_{i,i=1,...,n} & = \xi + \lambda           \sinh \left( \dfrac{\phi_i-\gamma }{\delta} \right) \\
	S_L : \eta_{i,i=1,...,n} & = \xi + \lambda            \exp \left( \dfrac{\phi_i-\gamma }{\delta} \right) \\
	S_B : \eta_{i,i=1,...,n} & = \xi + \lambda \left[ 1 + \exp \left( \dfrac{\phi_i-\gamma }{\delta} \right) \right]^{-1}
\end{align*}

The kind of transformation used (Unbounded, Log-normal or Bounded) and the \(\xi, \lambda, \gamma, \delta\) parameters depend on the prescribed moments.
Hill et al. **[16]** have provided an algorithm which automatically choose the right transformation with the associated parameters.

###Step 5
*\(z\) is obtained by digital filtering of \(\eta\)*

The digital filter \(H\) is used to average \(h\), and finally \(z = \text{FFT}^{-1}(Z=HA)\).
At the end of the process, \(z\) is supposed to exhibit the right autocorrelation function *acf* and the right four statistics, \(\mu\), \(\sigma\), \(Sk\) and \(Ku\).














