Q2D3
----

This is a prototype that illustrates how to use the QIIME 2 SDK to build an interface. This is based on the earlier [q2d2](https://github.com/gregcaporaso/q2d2) prototype.

As this is a prototype, you shouldn't rely on any of the code in this subpackage. This code will be removed in favor of production-ready interfaces in the near future.

Using this interface
--------------------

Run these install commands before you attempt to use this interface for the first time:

```bash
conda install notebook pyyaml
pip install https://github.com/rossant/ipymd/archive/master.zip
```

Then, to launch the server, change to this directory (``qiime2/q2d3``) and run the following command:

```bash
qiime-q2d3 serve
```
