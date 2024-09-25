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


Fractures have a great impact on oil exploration as they modify fluid flow properties within reservoir rocks, creating an interconnected network. Often, the oil reservoirs are difficult to assess, and remote sensing techniques offer too sparse data or have a low-spatial resolution. Otherwise, outcrops allow fracture characterization directly in the field or by digital representations in 2D and 3D. Fractures networks in this space can be used as source data providing fracture information that is generally scaled to reservoir modeling. This scaling is done by replicating the appropriate distribution and statistical attributes where length data is usually related to fractal propagation and power-law distribution parameters. In this sense, attribute estimators like the Maximum Likelihood Estimator (MLE) have been widely used presenting great approximation to actual parameters. However, our tested initial hypothesis found that optimum values for power-law distributions can be found arbitrarily despite the effectiveness of the MLE. In this sense, artificial intelligence and machine learning algorithms can be exploited for this kind of approximation. Our work proposes the use of an optimization algorithm based on Stochastic Gradient Descent (SGD) with momentum to obtain best-fitting parameters for power-law distributions. The proposed method is evaluated, firstly, with synthetic data and several metrics, with synthetic data against the traditional MLE, and at last, using empirical data obtained from fracture characterization in the Digital Outcrop Model (DOM) of a reservoir analogue outcrop. The results showed that the SGD method was able to provide better distribution fitting than the MLE and two other methods that also estimate the lower cutoff value. Albeit, slower than numerical methods, we tried to accelerate the algorithm execution by exploiting multiprocessing techniques. In conclusion, the optimization method is a valuable resource as a power-law attribute estimator.

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
