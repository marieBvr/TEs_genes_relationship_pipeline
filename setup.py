try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(
    name='te_abricot',
    version='1.0.0',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    py_modules=[
        'Create_Data',
        'Create_Data_LTR_multiprocessing',
        'Create_Data_LTR',
    ],
)
