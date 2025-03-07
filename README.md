<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
<!-- [![LinkedIn][linkedin-shield]][linkedin-url] -->



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/Ademartio/BPOTF/tree/OBPOTF">
    <img src=".github/imgs/BPOTF_logo.jpg" alt="Logo" width="160" height="160">
  </a>

  <h3 align="center">BPOTF Decoder</h3>

  <p align="center">
    Belief Propagation Ordered Tanner Forest decoder
    <br />
    <a href="https://github.com/Ademartio/BPOTF/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/Ademartio/BPOTF/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
    <li><a href="#attribution">Attribution</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This project implements a new decodig method for quantum low density parity check codes which we have named *Belief Propagation Ordered Tanner Forest* (BPOTF). The method is based on an slighlty modified Kruskal's algorithm to find the spanning forest associated to the columns of a PCM with highest a posteriori probability coming from running belief propagation and it uses the Disjoint-Set data structure for the implementation.

<!-- Add more info here... -->

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

This section explains which are the prerequites, how to compile and use the Python module created from this project.

### Prerequisites

The project has the following dependecies, which some of them are handled by pip/CMake installation processes and do not need to be accounted for. The dependencies are:

- C++20
- Pybind11 (automatic download)
- Cmake (pip installation supports it)
- [LDPC](https://github.com/quantumgizmos/ldpc.git) (part of its source code is embedded in the module)
- [BeliefMatching](https://github.com/oscarhiggott/BeliefMatching) (automatic handling with pip)
- [Stim](https://github.com/quantumlib/Stim) (automatic handling with pip)
- [SciPy](https://github.com/scipy/scipy) (automatic handling with pip)

### Installation

Below are the explanations on how to compile the python module to use the BPOTF decoder. It is highly recommended to use the pip option, while for developing use pip's editable mode and then the CMake option.

NOTE: the installation instructions are presented for linux systems, but some may also work for other kind of OS.

#### CMake 

The steps to compile the python importable BPOTF module using this method are explained below.

1. Clone the git repository and navigate to the directory.
   ```sh
   git clone https://github.com/Ademartio/BPOTF.git
   cd BPOTF
   ```
2. Create a new build folder and navigate to it.
   ```sh
   mkdir build
   cd build
   ```
3. Configure the Cmake project and build the python importable module.
   ```sh
   cmake ..
   make
   ```

When using this method alone, to be able to use the module it must be copied to the virtual environment or import the location of the module directly from the python script.

#### Pip

To compile the python importable module using pip, the only requisite is to have a C++20 compatible compiler and a pip 10+. Currently, Debian based systems and Windows 10 systems with MSVC 17 have been tested with this option. To compile it using this method, the steps indicated below can be followed:

1. Clone the git repository and navigate to the directory.
   ```sh
   git clone https://github.com/Ademartio/BPOTF.git
   cd BPOTF
   ```
2. (Optional) It is recommended to create a virtual environment to install python packages local to the project.
   ```sh
   python3 -m venv .venv
   source .venv/bin/activate # linux
   .venv\Scripts\activate # Windows
   ```
3. Execute the following command:
   ```sh
   pip install . # Normal install
   pip install -e . # Editable mode
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Simple usage of the module, import it just like another python module and use the offered OBPOTF object to initialize and decode the error. A minimal example of its initialization could be the following:

```Python
# Assuming the module is placed in the compilation output folder 'module'
from module import BPOTF
import numpy as np

_bpotf = BPOTF.OBPOTF(pcm.astype(np.uint8), p)

# ...
# Generate errors/syndrome or whatever
# ...

_bpotf.decode(syndrome.astype(np.int32))
```

<!-- _For more examples, please refer to the [Documentation](https://example.com)_ -->

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [x] Add a License to the repo.
- [x] Add CSC support.
- [x] Add Pip supported installation method.
- [ ] Add function to calculate an sparsified detector error model and transfer matrix.
- [x] Add function to map soft information from detector error model to sparsified detector matrix.
- [x] Add BP+BP decoding protocol for circuit-level noise.
- [ ] Add OTF function against LLRs (using logarithmic probabilities)
- [ ] Try to optimize BP decoder

See the [open issues](https://github.com/Ademartio/BPOTF/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

* _Antonio de Martì_ - [@ton_demarti](https://x.com/ton_demarti) - toni.demarti@gmail.com
* _Josu Etxezarreta_ - [@katutxakur](https://x.com/katutxakur) - jetxezarreta@unav.es
* _Imanol Etxezarreta_ [@ima96.bsky.social](https://bsky.app/profile/ima96.bsky.social) - ietxezarretam@gmail.com
* _Joschka Roffe_ - [@quantumgizmos](https://x.com/quantumgizmos) - joschka@roffe.eu


Project Link: [https://github.com/Ademartio/BPOTF](https://github.com/Ademartio/BPOTF)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Thanks to the following amazing projects and webs for the help, tools and information! Do not forget to visit and star/like their work also!

* [Pybind11](https://github.com/pybind/pybind11/tree/master)
* [LDPC python library](https://github.com/quantumgizmos/ldpc.git) - Joschka Roffe
* [Beliefmatching package](https://github.com/oscarhiggott/BeliefMatching) - Oscar Higgott
* [Stim](https://github.com/quantumlib/Stim)
* [SciPy](https://github.com/scipy/scipy)
* [Best-README-Template](https://github.com/othneildrew/Best-README-Template) - Othneil Drew

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ATTRIBUTION -->
## Attribution
When using the OTF post-processing decoding algorithm or the two stage BP decoder please cite our paper:
```
@article{bpotf_2024,
    author = "{deMarti iOlius}, Antonio and {Etxezarreta Martinez}, Imanol and Roffe, Joschka and {Etxezarreta Martinez}, Josu",
    title = "{An almost-linear time decoding algorithm for quantum LDPC codes under circuit-level noise}",
    journal = {arXiv},
    pages = {2409.01440},
    archivePrefix = "arXiv",
    primaryClass = "quant-ph",
    month = sep,
    year = {2024},
    url ={https://arxiv.org/abs/2409.01440}
}
```




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/Ademartio/BPOTF.svg?style=for-the-badge
[contributors-url]: https://github.com/Ademartio/BPOTF/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/Ademartio/BPOTF.svg?style=for-the-badge
[forks-url]: https://github.com/Ademartio/BPOTF/network/members
[stars-shield]: https://img.shields.io/github/stars/Ademartio/BPOTF.svg?style=for-the-badge
[stars-url]: https://github.com/Ademartio/BPOTF/stargazers
[issues-shield]: https://img.shields.io/github/issues/Ademartio/BPOTF.svg?style=for-the-badge
[issues-url]: https://github.com/Ademartio/BPOTF/issues
[license-shield]: https://img.shields.io/github/license/Ademartio/BPOTF.svg?style=for-the-badge
[license-url]: https://github.com/Ademartio/BPOTF/LICENSE
<!-- [linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew -->
