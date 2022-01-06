from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='gauchian',
      version='1.0',
      description='WGS-based GBA variant caller',
      long_description=readme(),
      classifiers=[],
      keywords='GBA',
      url='https://github.com/illumina/Gauchian',
      author='Xiao Chen',
      author_email='xchen2@illumina.com',
      license='GPLv3',
      packages=['gauchian', 'gauchian.caller', 'gauchian.depth_calling'],
      package_data={'gauchian': ['data/*']},
      install_requires=['pysam', 'numpy', 'scipy', 'statsmodels'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      entry_points={
          'console_scripts': [
              'gauchian=gauchian.gauchian:run'
              ]
      },
      include_package_data=True,
      zip_safe=False)
