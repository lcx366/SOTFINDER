from setuptools import setup,find_packages 

setup(
    name='sotfinder',
    version='0.0.1',
    description='SOTFINDER is an astrophotography analysis tool that expertly identifies and tracks space object trajectories in multi-frame images, while offering photometric analysis and dynamic visualization against a stellar backdrop.',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/SOTFINDER',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['Astrometric Analysis','Trajectory Tracking','Photometric Measurement','Space Object Detection'],
    python_requires = '>=3.10',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        'scipy',
        'matplotlib',
        'h5py',
        'starextractor',
        'scikit-learn',
        'photutils',
        'scikit-image'
        ],
)
