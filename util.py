# utilities for the web app

import numpy as np
import re
from lib_temper import *


def cents(x, prec=3):
	return '{1:.{0}f}'.format(prec, 1200 * x)


def parse_subgroup(s):
	s = [Fraction(i) for i in re.split("[\\.,; ]+", s)]

	if len(s) == 1:
		expanded = p_limit(s[0])
		return np.eye(len(expanded), dtype=np.int64), expanded
	else:
		s_basis, expanded = get_subgroup_basis(s)
		s = get_subgroup(s_basis, expanded)
		# s = expanded
		return s_basis, expanded


def parse_edos(s):
	s = [int(i) for i in re.split("[\\.,; &]+", s)]
	return s


ratioPattern = '(\d+)[/:](\d+)'
vectorPattern = '[\[(<]\s*(-?\d+(?:[,\s]+-?\d+)*)\s*[\])>]'


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

	commas = np.hstack(commas)

	M_expanded = hnf(cokernel(commas))

	commas_2 = solve_diophantine(basis, commas)

	M = hnf(cokernel(commas_2))

	assert np.allclose(M_expanded @ basis @ commas_2, 0), "comma not in basis"

	return (M, basis, M_expanded, s_expanded)


def from_edos(args):
	basis, s_expanded = parse_subgroup(args["subgroup"])
	edo_list = parse_edos(args["edos"])

	edos = []
	for e in edo_list:
		edos.append(patent_map(e, get_subgroup(basis, s_expanded)))

	M = hnf(np.vstack(edos), remove_zeros=True)

	# remove contorsion
	M = defactored_hnf(M)

	# find expansion with the same kernel
	M_expanded = hnf(cokernel(basis @ kernel(M)))

	return (M, basis, M_expanded, s_expanded)


def info(temp, options):
	T = temp[0]
	basis = temp[1]
	T_expanded = temp[2]
	s_expanded = temp[3]

	## we can find the map as:
	# T_expanded @ basis (removing zero rows)
	# so it is possibly redundant
	# however, it might be defactored!
	## example: 2.9.5.7/6 + 13&31
	# print(T)
	# print(T_expanded @ basis, flush = True)
	# print(basis, flush = True)
	# print(defactored_hnf(basis), flush = True)
	# print(T_expanded)
	# print("=========")
	# T = hnf(T_expanded @ basis, remove_zeros = True)

	s = get_subgroup(basis, s_expanded)

	res = dict()

	res["rank"] = T.shape[0]
	res["dim"] = T.shape[1]

	res["subgroup"] = ".".join(map(str, s))

	commas = LLL(kernel(T))

	comma_str = []
	for c in commas.T:
		comma_str.append(str(ratio(make_positive(c, s), s)))

	res["commas"] = ", ".join(comma_str)

	edolist = find_edos(T, s)
	print(edolist, flush = True)
	if edolist is not None and len(edolist) > 1:
		joins = find_join(T, s, edolist)

		show_list = [m[0][0][0] for m in edolist]

		res["edos"] = ', '.join(map(str, show_list))

		if joins is not None:
			res["edo joins"] = ' & '.join(map(str, joins))


		# res["edomapping"] = format_matrix(Tm)

	gens = preimage(T)

	if options["reduce"]:
		# eq = log_interval(gens[0], s)
		o = T[0, 0]
		genoct = np.zeros_like(gens[:, 0])
		genoct[0] = 1
		# reduce by octave
		for i in range(1, T.shape[0]):
			# make positive first
			if log_interval(gens[:, i], s) < 0:
				T[i, :] = -T[i, :]
				gens[:, i] = -gens[:, i]

			red = int(np.floor(log_interval(gens[:, i], s)))
			gens[:, i] -= red * genoct
			T[0, :] += o * red * T[i, :]

		# should be the same
		# gens = preimage(T)

	else:
		# make positive
		for i in range(T.shape[0]):
			if log_interval(gens[:, i], s) < 0:
				T[i, :] = -T[i, :]
				gens[:, i] = -gens[:, i]

	res["mapping"] = format_matrix(T)

	gens_print = [ratio(g, s) for g in gens.T]

	res["preimage"] = list(map(str, gens_print))

	weight = "unweighted"
	if options["tenney"]:
		weight = "tenney"

	te_tun, te_err = lstsq((T_expanded, s_expanded), weight)
	cte_tun, cte_err = cte((T_expanded, s_expanded), weight)

	te_tun = (te_tun.T @ T_expanded @ basis) @ gens
	cte_tun = (cte_tun.T @ T_expanded @ basis) @ gens

	res["TE-tuning"] = list(map(cents, te_tun.flatten()))

	res["CTE-tuning"] = list(map(cents, cte_tun.flatten()))
	res["TE-errors"] = ", ".join(map(cents, te_err))
	res["CTE-errors"] = ", ".join(map(cents, cte_err))

	# badness, complexity, error = temp_measures((T_expanded, s_expanded))
	# res["badness"] = '{:.3e}'.format(badness)
	# res["complexity"] = '{1:.{0}f}'.format(3, complexity)
	# res["error"] = '{1:.{0}f}'.format(3, 1200*error)

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