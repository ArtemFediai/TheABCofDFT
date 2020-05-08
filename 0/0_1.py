#! bin/bash

import os

a = os.path.dirname(os.path.abspath(__file__))
print('\nHello!')
print('I am the first script you are running at int-nano.')
print('My goal is to return the path to myself.')
print('This is it:', a)
print('Bye!!!\n')
