#import important library

import pandas as pd
import numpy as np
from rpy2.robjects import pandas2ri,r,Formula
from rpy2.robjects.packages import importr

# Import DESeq2
deseq2 = importr('DESeq2')

pandas2ri.activate()
