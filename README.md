# TorqTwin
This is the repository for uploading the codes for the following publication:<br>
TorqTwin—An open-source reference multibody modeling framework for wind turbine structural dynamics<br>
Feng Guo<sup>a</sup>,Zhen Gao<sup>a</sup>, David Schlipf<sup>b,c</sup> <br>
a. State Key Laboratory of Ocean Engineering, Shanghai Jiao Tong University, 800 Dongchuan Road, Shanghai, China<br>
b. Wind Energy Technology Institute, Flensburg University of Applied Sciences, Kanzleistraße 91-93, 24943 Flensburg, Germany<br>
c. sowento GmbH, Hessenlauweg 14, 70569 Stuttgart, Germany<br>
https://doi.org/10.1016/j.renene.2024.121268<br>

# Tutorial
To reproduce the simulation results presented in the paper, please do the following steps:

1. Create a virtual environment by:
   <pre>conda create -c conda-forge -n TorqTwin python=3.10 (recommanded)</pre>
2. Activate the virtual environment
   <pre>conda activate TorqTwin </pre>
3. Install required packages
   <pre>conda install -c conda-forge ipympl ipython jupyter notebook matplotlib numpy pythreejs "scikits.odes" scipy "sympy>=1.12"</pre>
4. Open Jupyter Notebook
   <pre>jupyter notebook</pre>
5. Run examples
   Use jupyter notebook browser to open an example file, e.g. "ElastoDyn_fixbottom/TorqTwin_Onshore.ipynb", and run the cells.
   


