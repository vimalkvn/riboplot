# -*- coding: utf-8 -*-
import os


class Config(object):
    # get the running directory of this file, move one level up to get the
    # application directory
    APP_DIR = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    PKG_DATA_DIR = os.path.join(APP_DIR, 'riboplot', 'data')


class TestingConfig(Config):
    """Testing configuration"""
    TEST_DATA_DIR = os.path.join(Config.APP_DIR, 'tests/data')


class ProductionConfig(Config):
    """Production configuration"""
    # Additional variables can be listed here
    pass
