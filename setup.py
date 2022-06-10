from setuptools import setup, find_packages

setup(name='PVSdock',
        version='0.1',
        packages=['pvsdock'],
#        packages=find_packages(),
        url='https://github.com/gicsaw/PVSdock',
        license='MIT LICENSE',
        author='Seung Hwan Hong',
        author_email='gicsaw0@gmail.com',
        description='',
        scripts=['bin/master_dock.py',
                 'bin/pydock_run.py',
                 'bin/sub_dock.py',
                 'bin/vhts_check_restart.py'
                ]
)


