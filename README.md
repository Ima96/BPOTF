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
  <a href="https://github.com/Ima96/BPOTF">
    <img src=".github/imgs/BPOTF_logo.jpg" alt="Logo" width="160" height="160">
  </a>

  <h3 align="center">BPOTF Decoder</h3>

  <p align="center">
    Belief Propagation Ordered Tanner Forest decoder
    <br />
    <a href="https://github.com/Ima96/BPOTF/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/Ima96/BPOTF/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
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
    <li><a href="#package-content">Package Content</a></li>
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
- [Stim](https://github.com/quantumlib/Stim) (automatic handling with pip)
- [SciPy](https://github.com/scipy/scipy) (automatic handling with pip)
- [Sinter](https://pypi.org/project/sinter/) (automatic handling with pip)

### Installation
Below are the explanations on how to compile the python module to use the BPOTF decoder. It is highly recommended to use the pip option, while for developing use pip's editable mode and then the CMake option.

*NOTE: the installation instructions are presented for linux systems, but some may also work for other kind of OS.*

<details>
<summary><b>CMake</b></summary>

This method is thought for developers that want to improve or debug the code, so that they do not need to execute the whole pip process for every little change. Thus, install first the package using pip in editable mode, and then use this workflow.

The steps to compile the python importable BPOTF module using this method are explained below.

1. Clone the git repository and navigate to the directory.
   ```sh
   git clone https://github.com/Ima96/BPOTF.git
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

</details>

<details>
<summary><b>Pip</b></summary>

To compile the python importable module using pip, the only requisite is to have a C++20 compatible compiler and a pip 10+. Currently, Debian based systems and Windows 10 systems with MSVC 17 have been tested with this option. To compile it using this method, the steps indicated below can be followed:

1. Clone the git repository and navigate to the directory.
   ```sh
   git clone https://github.com/Ima96/BPOTF.git
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

</details>

<details>
<summary><b>PyPI (Coming soon)</b></summary>

As any other package:
```sh
pip install BPOTF
```

</details>

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- PACKAGE EXPLANATION -->
## Package Content
The Python package exposes the following objects from the `BPOTF` module: the `OBPOTF` decoder class, the `DemData` helper struct, the `SinterBpOtfDecoder` sinter adapter class, and the `NoiseType` enumeration that controls the decoding mode.

---

### `NoiseType` enumeration
Selects the noise model and decoding pipeline used internally:

| Value | Integer | Description |
|-------|---------|-------------|
| `NoiseType.E_CC` | 0 | **Code Capacity** (default). Single-stage BP followed by the OTF spanning-forest post-processing. Use this for standard PCM-based simulations without circuit-level noise. |
| `NoiseType.E_PHEN` | 1 | **Phenomenological**. Reserved for phenomenological noise models. Not yet supported. |
| `NoiseType.E_CLN` | 2 | **Circuit-Level Noise**. Enables a two-stage BP+BP+OTF pipeline. Requires a populated `DemData` object to supply the sparsified detector error model matrices. |

---

### `DemData` struct
A helper data container used when constructing an `OBPOTF` object in **Circuit-Level Noise** (`E_CLN`) mode. All fields accept NumPy arrays or SciPy sparse matrices and default to `None`.

| Attribute | Description |
|-----------|-------------|
| `obs_matrix` | Observable matrix derived from the DEM (CSR format). |
| `phen_obs_matrix` | Sparsified observable matrix (CSR format). |
| `phen_check_matrix` | Sparsified parity-check matrix (CSC format). Used as the first-stage PCM in the two-stage BP decoder. |
| `transfer_matrix` | Transfer matrix that maps soft information from the detector-error-model space to the sparsified detector matrix (CSR format). |
| `priors` | 1-D array of prior error probabilities for the detector error model columns. |

---

### `OBPOTF` class
The main decoder object. The constructor selects the internal decoding callback based on `noise_type` and the data provided.

#### Constructor
```python
BPOTF.OBPOTF(
    pcm,                          # Parity-check matrix (np.uint8 NumPy array or SciPy CSC matrix)
    p,                            # Physical error probability (float)
    noise_type = NoiseType.E_CC,  # Noise model / decoding mode
    po_otf_csc_mat = None,        # Optional: separate matrix on which to apply the OTF step
    po_ext_bp_iters = None,       # Optional: array of three ints [pcm_bp, phen_bp, otf_bp] to override default BP iteration counts
    ps_ext_dem_data = None,       # Optional: DemData object (required for E_CLN two-stage mode)
    decimation = 1e-9             # Probability assigned to columns not selected by OTF
)
```

**Decoding modes selected by `noise_type`:**
- **`E_CC` — Code-Capacity mode** (default)  
  A single BP pass is run over `pcm`, and the resulting posterior probabilities drive the OTF spanning-forest algorithm. Best suited for depolarising or i.i.d. noise without measurement errors.

- **`E_CLN` without `DemData`**  
  BP is run over the circuit-level PCM and the OTF post-processing is applied directly to it. No phenomenological pre-processing is performed.

- **`E_CLN` with `DemData` (two-stage BP+BP+OTF)**  
  Activates the full two-stage pipeline. A first BP pass is run over the phenomenological PCM (`phen_check_matrix`) to compute soft information, which is then propagated through the `transfer_matrix` to produce priors for a second BP pass over the main circuit-level PCM. The OTF step is finally applied to the second-stage output. This mode yields the best performance under circuit-level noise.

#### Methods
| Method | Description |
|--------|-------------|
| `decode(syndrome)` | Runs the configured decoding pipeline on the given syndrome vector (`np.uint8` array). Returns the recovered error as a `np.uint8` array. |
| `has_converged()` | Returns `True` if the most recent BP pass converged before reaching its iteration limit. |
| `print_object()` | Prints the object's internal state. Intended for debugging. |

---

### BP iteration defaults
When `po_ext_bp_iters` is not provided, the following defaults are used internally:

| Stage | Default iterations |
|-------|--------------------|
| First-stage PCM BP (`pcm_bp_iters`) | 30 |
| Phenomenological BP (`phen_bp_iters`) | 100 |
| Post-OTF BP (`otf_bp_iters`) | 100 |

Pass a 3-element integer array to `po_ext_bp_iters` to override all three values.

---

### `SinterBpOtfDecoder` class
A [Sinter](https://pypi.org/project/sinter/) adapter that wraps `OBPOTF` so it can be plugged directly into `sinter.collect()` for benchmarking under circuit-level noise. It implements `sinter.Decoder` and always uses the full **BP+BP+OTF** pipeline (`E_CLN` with `DemData`).

> **Note:** The `OBPOTF` object is built *lazily* on the first call to `decode_via_files` / `decode`, because Sinter serialises decoder objects across worker processes and the underlying C++ objects are not picklable.

#### Constructor
```python
from BPOTF import SinterBpOtfDecoder

decoder = SinterBpOtfDecoder(
    dem,                    # stim.DetectorErrorModel — the only mandatory argument
    p,                      # Physical error probability (float)
    pcm              = None,  # Override PCM (SciPy sparse). Derived from DEM if None.
    obs_matrix       = None,  # Override observable matrix. Derived from DEM if None.
    transfer_matrix  = None,  # Transfer matrix (enables three-stage path when provided).
    phen_check_matrix= None,  # Sparsified PCM. Derived from DEM edge_check_matrix if None.
    phen_obs_matrix  = None,  # Sparsified observable matrix. Derived from DEM if None.
    otf_matrix       = None,  # Matrix on which the OTF step is applied. Defaults to sparsified PCM.
    bp_iters         = None,  # [stage1, stage2, stage3] int list. Uses OBPOTF defaults if None.
    decimation       = 1e-9   # OTF decimation parameter.
)
```

When `transfer_matrix` is `None` the decoder still builds a valid `DemData` object but the two-stage soft-information propagation step is skipped internally by `OBPOTF`.

All matrices that are not explicitly provided are derived automatically from the `stim.DetectorErrorModel` via `ldpc.ckt_noise.dem_matrices.detector_error_model_to_check_matrices`. The priors are always taken from the DEM regardless.

#### Methods
| Method | Description |
|--------|-------------|
| `decode(syndrome)` | Decodes a single syndrome vector. Lazily builds the `OBPOTF` object on first call. Returns a `np.ndarray` of observable predictions. |
| `decode_via_files(...)` | Sinter worker entry-point. Reads shot data from `dets_b8_in_path`, runs `decode` for every shot, and writes observable predictions to `obs_predictions_b8_out_path`. |

#### Example usage with Sinter
There is an Jupyter Notebook showcasing an example usage of the sinter decoder. It can be found in `examples/sinter_example.ipynb`.


<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/Ima96/BPOTF/issues) for a full list of proposed features (and known issues).

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
* _Antonio de Martí_ - [@ton_demarti](https://x.com/ton_demarti) - toni.demarti@gmail.com
* _Josu Etxezarreta_ - [@katutxakur](https://x.com/katutxakur) - jetxezarreta@unav.es
* _Imanol Etxezarreta_ - ietxezarretam@gmail.com
* _Joschka Roffe_ - [@quantumgizmos](https://x.com/quantumgizmos) - joschka@roffe.eu


Project Link: [https://github.com/Ima96/BPOTF](https://github.com/Ima96/BPOTF)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
Thanks to the following amazing projects and webs for the help, tools and information! Do not forget to visit and star/like their work also!

* [Pybind11](https://github.com/pybind/pybind11/tree/master)
* [LDPC python library](https://github.com/quantumgizmos/ldpc.git) - Joschka Roffe
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
[contributors-shield]: https://img.shields.io/github/contributors/Ima96/BPOTF.svg?style=for-the-badge
[contributors-url]: https://github.com/Ima96/BPOTF/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/Ima96/BPOTF.svg?style=for-the-badge
[forks-url]: https://github.com/Ima96/BPOTF/network/members
[stars-shield]: https://img.shields.io/github/stars/Ima96/BPOTF.svg?style=for-the-badge
[stars-url]: https://github.com/Ima96/BPOTF/stargazers
[issues-shield]: https://img.shields.io/github/issues/Ima96/BPOTF.svg?style=for-the-badge
[issues-url]: https://github.com/Ima96/BPOTF/issues
[license-shield]: https://img.shields.io/github/license/Ima96/BPOTF.svg?style=for-the-badge
[license-url]: https://github.com/Ima96/BPOTF/LICENSE
<!-- [linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew -->
