# SGD_power_law_estimator
Power-law distribution attribute estimator based on SGD. This is the code and data companion for the following paper:

**[Paper] Stochastic Gradient Descent optimization to estimate the power-law fractal index in fracture networks**


Authors:
 [Graciela dos Reis Racolte<sup>1</sup>](https://www.researchgate.net/profile/Graciela-Racolte),
 [Ademir Marques Junior<sup>1</sup>](https://www.researchgate.net/profile/Ademir_Junior),
 [Eniuce Menezes<sup>2</sup>](https://www.researchgate.net/profile/Eniuce-Souza),
 [Leonardo Scalco<sup>1</sup>](https://www.researchgate.net/profile/Leonardo-Scalco),
 [Delano Ibanez<sup>3</sup>](),	
 [Maurício Roberto Veronez<sup>1</sup>](https://www.researchgate.net/profile/Mauricio_Veronez),
  [Luiz Gonzaga Junior<sup>1</sup>](https://www.researchgate.net/profile/Luiz_Gonzaga_da_Silveira_Jr)
 
[Vizlab | X-Reality and GeoInformatics Lab<sup>1</sup>](http://vizlab.unisinos.br/), 
[State University of Maringá, Paraná, Brazil<sup>2</sup>](http://www.cpr.uem.br/international/index.php/en/),
[CENPES-Petrobras<sup>3</sup>](https://petrobras.com.br/en/our-activities/technology-innovation/)  
Submitted to Elsevier's Computer & Geoscience


Fractures greatly impact hydrocarbon exploration as they modify fluid flow properties within reservoir rocks, creating an interconnected network. The hydrocarbon reservoirs are often difficult to assess, and the methods employed in acquiring information from these locations offer too sparse data or have a low spatial resolution. Otherwise, outcrops allow fracture characterization directly in the field or using 2D and 3D digital representations of outcrops. These fracture networks, usually related to fractal propagation and power-law distribution parameters, can be used as data sources providing useful information when properly adjusted to the reservoir simulation scale. In this sense, attribute estimators, like the Maximum Likelihood Estimator (MLE) and algorithms using MLE, have been widely used for their robustness when compared to linear regression estimators. However, due to the challenges in the power-law characterization, such as the large fluctuations that occur in the tail of the distribution, non-optimum values can be obtained despite the effectiveness of the MLE. Our work proposes the use of an optimization algorithm based on Stochastic Gradient Descent (SGD) with momentum to obtain best-fitting parameters for power-law distributions. The proposed method was first evaluated with synthetic data and several goodness-of-fitness metrics and later using empirical data obtained from fracture characterization in the Digital Outcrop Model (DOM) of a reservoir analogue outcrop. Stochastic DFN sampling based on empirical data was also used to simulate censoring effects. The results showed that the SGD method provided better distribution fitting than other methods based on the MLE when using empirical data while presenting reduced bias when using synthetic data. The estimation of power-law parameters in stochastic DFN data also presented the best-fitting results when applying the proposed method. In conclusion, the proposed optimization method proved a valuable alternative to estimate power-law distributions.

## Installation

This Python script requires a Python 3 environment and the following installed libraries as seen in the file requirements.txt.

## Usage

Use the link below to open the script on Google Colab.

[![Google Colab](https://badgen.net/badge/Launch/on%20Google%20Colab/blue?icon=terminal)](https://colab.research.google.com/drive/1DPhNua2DDCjqB8ZMEFNuIGP6DWMOcR8Y?usp=sharing)

Fracture interpretation used in the paper is provided in this diretory and need to be placed in the same directory of the code if run locally, or placed in the Colab environment or Google drive.

## Credits	
This work is credited to the [Vizlab | X-Reality and GeoInformatics Lab](http://vizlab.unisinos.br/) and the following developers:	[Ademir Marques Junior](https://www.researchgate.net/profile/Ademir_Junior) and [Graciela dos Reis Racolte](https://www.researchgate.net/profile/Graciela-Racolte).

## License

    MIT Licence (https://mit-license.org/)

## How to cite

```bash
@article{racolte2024stochastic,
  title={Stochastic Gradient Descent optimization to estimate the power-law fractal index in fracture networks},
  author={Racolte, Graciela and Marques Jr, Ademir and Menezes, Eniuce and Scalco, Leonardo and Ibanez, Delano Menecucci and Veronez, Mauricio Roberto and Gonzaga Jr, Luiz},
  journal={Computers \& Geosciences},
  volume={191},
  pages={105677},
  year={2024},
  publisher={Elsevier}
}
```
