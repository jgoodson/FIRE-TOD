# Using python on TOD-compute

Although you are free to install Python on your own machine and work locally, I have prepared a Python installation with a large variety of useful packages pre-installed on the TOD-compute server. To access this, you will need to run a single command:

```bash
/opt/anaconda/bin/conda init
```

This will update your user account to access the Anaconda install of Python and set your default Python to the Anaconda install of Python 3.7. After you run this you will need to log out and back in. Once you do so, if you run 

```bash
which python
```

you should see the output `/opt/anaconda/bin/python`. With that, you should be set up to access our base install and you will not need to install package like numpy, scikit-learn, biopython, or more!