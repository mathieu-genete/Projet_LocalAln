#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:25:24 2023

@author: mathieu
"""

__title__ = 'AlignLib'
__version__ = '1.0'
__author__ = 'Mathieu GENETE'
__license__ = 'creative commons'
__copyright__ = 'Copyright 2024 Mathieu GENETE'


from .Alphabet import Alphabet
from .FastaSeq import FastaSeq
from .Seq import Seq
from .FilesParser import FilesParser
from .Cost import Cost
from .Alignment import Alignment
from .IUPAC import IUPAC
from .Output_screen import Output_screen

__all__ = [
        "Alphabet",
        "FastaSeq",
        "Seq",
        "FilesParser",
        "Cost",
        "Alignment",
        "IUPAC",
        "Output_screen"]