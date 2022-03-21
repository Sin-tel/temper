# utilities for the web app

import numpy as np
import re
from temper import *


def cents(x, prec=3):
	return '{1:.{0}f}'.format(prec, 1200 * x)


def tlist(l):
	if type(l) is list:
		return l
	else:
		return [l]


def parse_subgroup(s):
	s = [int(i) for i in re.split("[\\.,; ]+", s)]

	if len(s) == 1:
		return p_limit(s[0])
	else:
		return s


def parse_edos(s):
	s = [int(i) for i in re.split("[\\.,; &]+", s)]
	return s


ratioPattern = '(\d+)[/:](\d+)'
vectorPattern = '[[(<]\s*(-?\d+(?:[,\s]+-?\d+)*)\s*[])>]'


def parse_commas(c, s):
	commas = []
	for n, d in re.findall(ratioPattern, c):
		commas.append(factors((int(n), int(d)), s))
	for v in re.findall(vectorPattern, c):
		l = len(s)
		res = np.zeros((l, 1), dtype=np.int64)
		v = np.array(list(map(int, v.replace(',', ' ').split())))
		res[:v.shape[0], 0] = v[:l]
		commas.append(res)

	return commas


def format_matrix(matrix):
	"""Format a matrix using LaTeX syntax"""
	body_lines = [" & ".join(map(str, row)) for row in matrix]
	body = "\\\\\n".join(body_lines)
	body = "\\begin{bmatrix}" + body + "\\end{bmatrix}"
	return "\\( " + body + " \\)"


def from_commas(args):

	subgroup = parse_subgroup(args["subgroup"])

	commas = parse_commas(args["commas"], subgroup)

	res = hnf(cokernel(np.hstack(commas)))

	return (res, subgroup)


def from_edos(args):
	subgroup = parse_subgroup(args["subgroup"])
	edo_list = parse_edos(args["edos"])

	edos = []
	for e in tlist(edo_list):
		edos.append(patent_map(e, subgroup))

	res = hnf(np.vstack(edos), remove_zeros=True)

	# remove contorsion, for now
	res = defactored_hnf(res)

	return (res, subgroup)


def info(temp, options):
	T = temp[0]
	s = temp[1]

	res = dict()

	res["rank"] = T.shape[0]
	res["dim"] = T.shape[1]

	res["subgroup"] = ".".join(map(str, s))

	Tm = findMaps(T, s)

	res["edos"] = ' & '.join(map(str, list(Tm[:, 0])))

	# res["contorsion"] = factor_order(T)

	gens = preimage(T)

	if options["reduce"]:
		# eq = log_interval(gens[0], s)
		o = T[0, 0]
		genoct = np.zeros_like(gens[0])
		genoct[0] = 1
		# reduce by octave
		for i in range(1, T.shape[0]):
			# make positive first
			if log_interval(gens[i], s) < 0:
				T[i, :] = -T[i, :]
				gens[i] = -gens[i]

			print(log_interval(gens[i], s))
			red = int(np.floor(log_interval(gens[i], s)))
			gens[i] -= red * genoct
			T[0, :] += o * red * T[i, :]

		print(gens)
		# gens = preimage(T) # should be the same

	else:
		# make positive
		for i in range(T.shape[0]):
			if log_interval(gens[i], s) < 0:
				T[i, :] = -T[i, :]
				gens[i] = -gens[i]

	commas = LLL(kernel(T))

	comma_str = []
	for c in commas.T:
		comma_str.append(str(ratio(make_positive(c, s), s)))

	res["commas"] = ", ".join(comma_str)

	res["mapping"] = format_matrix(T)

	gens_print = [ratio(g, s) for g in gens]
	# res["preimage"] = ", ".join(map(str,gens_print))
	res["preimage"] = list(map(str, gens_print))

	weight = "unweighted"
	if options["tenney"]:
		weight = "tenney"

	te_tun, te_err = lstsq(temp, weight)
	cte_tun, cte_err = cte(temp, weight)
	res["te-tuning"] = list(map(cents, te_tun))

	res["cte-tuning"] = list(map(cents, cte_tun))
	res["te-error"] = ", ".join(map(cents, te_err))
	res["cte-error"] = ", ".join(map(cents, cte_err))

	return res
