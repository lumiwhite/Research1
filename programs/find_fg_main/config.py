import numpy as np
import random
import math
import sympy as sp
import itertools

from dataclasses import dataclass
import matplotlib.pyplot as plt
from collections import Counter
import sympy as sp
from sympy import gcd, gcdex, Matrix, list2numpy, eye
from itertools import product
from z3 import *

P = 384
J = 3
L = 12
l_h = L//2