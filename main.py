import numpy as np
import re
from temper import *

## todo
# mapping Gens

# print("Welcome to tuning calc v0.0.1 !")
# print("Type 'help()' for help.")

# def help():
#   print("To get info from joining edos, use 'from_edos(edos, subgroup)'.")
#   print("'edos' can be a single integer or a list in square brackets")
#   print("'subgroup' is a list of primes. If you put a single number it will take that prime limit.")
#   print("examples:")
#   print(">> from_edos([12,22], 7)\t\t(pajara)")
#   print(">> from_edos([5,17], [2,3,7])\t\t(archy)")
#   print("To temper out commas, use 'from_commas(commas, subgroup)'.")
#   print("'commas' is a single or a list of commas. They have to be in quotes!")
#   print("'subgroup', as above.")
#   print("examples:")
#   print('>> from_commas(["81/80","126/125"], 7)\t\t(septimal meantone)')
#   print('>> from_commas("225/224", 7)\t\t(This gives marvel)')


def cents(x):
  return '{0:.2f}'.format(1200*x)

def tlist(l):
  if type(l) is list:
    return l
  else:
    return [l]

def parse_subgroup(s):
  print("test")
  print(s)
  print(re.split("[\\.,; ]+",s))
  pattern = '-+'
  string = '2344------HELLO--WORLD'
  print(re.split(pattern, string))
  print("s")

  s = [int(i) for i in re.split("[\\.,; ]+",s)]

  if len(s) == 1:
    return p_limit(s[0])
  else:
    return s

def parse_edos(s):
  s = [int(i) for i in re.split("[\\.,; &]+",s)]
  return s

def from_commas(comma_list, s):

  subgroup = parse_subgroup(s)

  commas = []
  for l in tlist(comma_list): 
    commas.append(factors(l,subgroup))

  res = hnf(cokernel(np.hstack(commas)))

  return (res, subgroup)

def format_matrix(matrix):
  """Format a matrix using LaTeX syntax"""
  body_lines = [" & ".join(map(str, row)) for row in matrix]
  body = "\\\\\n".join(body_lines)
  body = "\\begin{bmatrix}" + body + "\\end{bmatrix}" 
  return "\\( " + body + " \\)" 

def from_edos(args):

  subgroup = parse_subgroup(args["subgroup"])
  edo_list = parse_edos(args["edos"])

  edos = []
  for e in tlist(edo_list): 
    edos.append(patent_map(e,subgroup))

  res = hnf(np.vstack(edos), remove_zeros = True)

  return (res, subgroup)



def info(temp):
  T = temp[0]
  s = temp[1]

  res = dict()

  res["subgroup"] = ".".join(map(str, s))

  Tm = findMaps(T, s)

  res["edos"] = ' & '.join(map(str, list(Tm[:,0])))

  res["contorsion"] = factor_order(T)

  commas = LLL(kernel(T))

  comma_str = []
  for c in commas.T:
    comma_str.append(str(ratio(c, s)))

  res["commas"] = ", ".join(comma_str)

  res["mapping"] = format_matrix(T)

  te_tun, te_err = lstsq(temp)
  cte_tun, cte_err = cte(temp)
  res["te-tuning"] = ", ".join(map(cents,te_tun))
  res["te-error"]  = ", ".join(map(cents,te_err))

  res["cte-tuning"] = ", ".join(map(cents,cte_tun))
  res["cte-error"]  = ", ".join(map(cents,cte_err))

  return res


def lstsq(temp, weight = "tenney"):
  M = temp[0]
  s = temp[1]

  j = np.log2(np.array(s))

  W = np.diag(1./j)
  if weight == 'unweighted':
    W = np.eye(len(s))

  j = np.atleast_2d(j)


  sol = np.linalg.lstsq((M@W).T, (j@W).T, rcond=None)[0]

  tun = (sol.T @ M)
  err = tun - j

  return sol.flatten(), err.flatten()

def cte(temp, weight = "tenney"):
  M = temp[0]
  s = temp[1]

  j = np.log2(np.array(s))

  W = np.diag(1./j)
  if weight == 'unweighted':
    W = np.eye(len(s))

  j = np.atleast_2d(j)

  A = (M@W).T
  b = (j@W).T

  r, d = M.shape

  V = np.zeros((d,1))
  V[0,0] = 1

  C = (M@V).T
  d = (j@V).T

  sol = np.linalg.solve(np.block([[A.T @ A, C.T],[C, 0]]), np.vstack([A.T @ b, d]))

  sol = sol[:r]

  tun = (sol.T @ M)
  err = tun - j

  return sol.flatten(), err.flatten()


# A = from_commas(["81/80","126/125"], 7)
# A = from_commas("81/80", 7)
# A = from_edos([5,19], 5)

# info(A)