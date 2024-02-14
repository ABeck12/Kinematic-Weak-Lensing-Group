**Fast Initial REdshift Fitting of cLuster galaxY (FIREFLY)**


To run, use your bash terminal and type:


.. code-block:: RST

    $ python main_firefly.py input_spec1d/ -all vvds_spiral


Here, ``input_spec1d/`` is the folder that you want FIREFLY to fit for redshifts, ``-all`` means you'd like to fit all files in your specified folder, ``vvds_spiral`` is the template that FIREFLY will use during the fitting.


To redo fitting for one spectrum, e.g., ``spec1d.m46.030.A2552.fits``, (and with your desired redshift, for example, z=0.6724), use 


.. code-block:: RST

    $ python main_firefly.py input_spec1d/ spec1d.m46.030.A2552.fits vvds_spiral 0.6724


See the updated documentation at `Google Colab <https://colab.research.google.com/drive/1s5pAIuA5Ou4Olkoos1lXTkWuoDD_Zf_d?usp=sharing>`_.
