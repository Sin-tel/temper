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
	s = [Fraction(i) for i in re.split("[\\.,; ]+", s)]

	if len(s) == 1:
		return p_limit(s[0])
	else:
		s_basis, expanded = get_subgroup_basis(s)
		s = get_subgroup(s_basis, expanded)
		# s = expanded
		return s_basis, expanded


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

	basis, s_expanded = parse_subgroup(args["subgroup"])

	commas = parse_commas(args["commas"], s_expanded)

	M_expanded = hnf(cokernel(np.hstack(commas)))

	# find commad expressed in subgroup basis
	# this is algorithm not exact!
	R = basis @ basis.T
	Rdet = np.linalg.det(R)
	Rinv = (np.linalg.inv(R) * np.linalg.det(R))
	Rinv = (np.round(Rinv).astype(np.int64))

	commas_2 = []
	for c in commas:
		sol = (Rinv @ basis @ c) / Rdet

		assert np.allclose(sol, np.round(sol)), "Comma not in subgroup"
		sol = np.round(sol).astype(np.int64)
		commas_2.append(sol)

	M = hnf(cokernel(np.hstack(commas_2)))

	assert np.allclose(M_expanded @ basis.T @ np.hstack(commas_2), 0)

	return (M, basis, M_expanded, s_expanded)


def from_edos(args):
	basis, s_expanded = parse_subgroup(args["subgroup"])
	edo_list = parse_edos(args["edos"])

	edos = []
	for e in tlist(edo_list):
		edos.append(patent_map(e, get_subgroup(basis, s_expanded)))

	M = hnf(np.vstack(edos), remove_zeros=True)

	# remove contorsion
	M = defactored_hnf(M)

	# find expansion from subgroup
	M_expanded = hnf(cokernel(basis.T @ kernel(M)))

	return (M, basis, M_expanded, s_expanded)


def info(temp, options):
	T = temp[0]
	basis = temp[1]
	T_expanded = temp[2]
	s_expanded = temp[3]

	s = get_subgroup(basis, s_expanded)

	res = dict()

	res["rank"] = T.shape[0]
	res["dim"] = T.shape[1]

	res["subgroup"] = ".".join(map(str, s))

	Tm = findMaps(T, s)

	res["edos"] = ' & '.join(map(str, list(Tm[:, 0])))

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

			red = int(np.floor(log_interval(gens[i], s)))
			gens[i] -= red * genoct
			T[0, :] += o * red * T[i, :]

		# should be the same
		# gens = preimage(T)

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
	# print(gens_print)

	res["preimage"] = list(map(str, gens_print))

	weight = "unweighted"
	if options["tenney"]:
		weight = "tenney"

	g_matrix = np.vstack(gens).T

	te_tun, te_err = lstsq((T_expanded, s_expanded), weight)
	cte_tun, cte_err = cte((T_expanded, s_expanded), weight)

	te_tun = (te_tun.T @ T_expanded @ basis.T) @ g_matrix
	cte_tun = (cte_tun.T @ T_expanded @ basis.T) @ g_matrix

	res["te-tuning"] = list(map(cents, te_tun.flatten()))

	res["cte-tuning"] = list(map(cents, cte_tun.flatten()))
	res["te-error"] = ", ".join(map(cents, te_err))
	res["cte-error"] = ", ".join(map(cents, cte_err))

	return res


#########################################
if __name__ == '__main__':
	args = dict()
	# args["subgroup"] = "2.3.5.7"
	args["subgroup"] = "2.5.9/7"
	# args["edos"] = "19,22"
	args["commas"] = "225/224"

	options = dict()

	options["tenney"] = True
	options["reduce"] = True

	# temp = from_edos(args)
	temp = from_commas(args)
	html_info = info(temp, options)