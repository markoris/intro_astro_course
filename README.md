# A Brief Introduction to Astrophysics

This repository contains code notebooks and student manuals covering introductory astrophysics material aimed at high-school level students. The content is designed in a modular fashion such that individual topics can be taught without dependence on other content. The instructor manual can be found [here](https://ccrgpages.rit.edu/~markoristic/intro_astro_materials/intro_astro_teaching_manual.pdf).

The code notebooks run locally on inexpensive (~$200) 14-inch HP Stream Laptops (Hardware Model ``HP HP Laptop 14s-cf2xxx") laptops running Ubuntu 22.04.2 LTS 64-bit. Hardware costs can, in theory, be further minimized to even more inexpensive Chromebooks if using the Google Colab cloud computing option (see below), although we have not personally tested this possibility. The instructions for setting up laptops in a similar fashion to ours can be found [here](https://ccrgpages.rit.edu/~markoristic/intro_astro_materials/intro_astro_laptop_setup.pdf).

# Code Notebook Usage

## Running on Google Colab (Preferred, First-Time-User Friendly)

We have created a [Google Colab notebook](https://colab.research.google.com/drive/1Gyz5FoMfnlCQS9s8Cu3zhnVXs8fZ_Kvw?usp=sharing) for easy interaction using Google Colab's cloud computing resources. This option only requires users to have a free Google account and should run smoothly on even the most computationally limited machines. 

After opening the notebook link above, simply click "Copy to Drive" in the top left to copy a version of the notebook to your own Google account. This copied version can be edited and saved, making it an excellent option for classroom settings where students can experiment with the code and save the edited notebooks.

## Running Locally (More Advanced)

The following guide assumes working knowledge of git repositories. Users are free to clone this repository to a local machine rather than using Google Colab's resources. This is a particularly good option if users anticipate adding their own content to the teaching material; in that case, we especially recommend creating a fork of this repository to have a version-controlled repository of your own. If you do add content of your own, please reach out to us and let us know! We would love to see your additions and incorporate them into our repository.

Once the repository is cloned locally, follow [these instructions](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) to set up a Python virtual environment.
Once finished with the [activate a virtual environment](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#activate-a-virtual-environment) step, use the requirements file to [install the required Python dependencies](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#using-a-requirements-file).

With the virtual environment setup and activated, the notebooks should run smoothly using Jupyter notebooks ([tutorial](https://www.codecademy.com/article/how-to-use-jupyter-notebooks)). To start, navigate to the intro_astro_course directory and run `jupyter notebook`. From there, you should be able to navigate to the respective sub-directories (e.g. instrumentation, observational, theoretical) and open the respective `.ipynb` to start using the code in the notebook.

# Workshop Website

You can find Kate's website about the course [here](https://kjtwagner.github.io/intro_astro_course/index.html), including [pictures](https://kjtwagner.github.io/intro_astro_course/gallery.html) from some of our previous workshops!

# Contact Us!

If you have any questions about the content, installation, usage, or contribution to any of the content presented here, please do not hesitate to reach us at [mr3364@rit.edu](mailto:mr3364@rit.edu) and [kjw4822@rit.edu](mailto:kjw4822@rit.edu)!
