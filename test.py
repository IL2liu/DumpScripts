#! /usr/bin/env python
# -*- coding: utf-8 -*-

# test.py

wanted = None

# caller.py
from test import wanted

'''
many prerequisites
'''
def imp_now(case):
    import test
    if case:
        test.wanted = test.case

pack = imp_now("limits_mm10")
print(pack)

