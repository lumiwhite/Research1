import math
import itertools
import numpy as np
import random
from dataclasses import dataclass
import matplotlib.pyplot as plt
from collections import Counter
import sympy as sp
from sympy import gcd, gcdex, Matrix, list2numpy, eye
from itertools import product
from z3 import *

L=12
l_h = L // 2
J = 3
P=768