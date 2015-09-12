import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def add_chord(ax, R, angle_a, angle_b, color, width):
	# cartesian coordinates of endpoint a
	x_a = R*np.cos(angle_a); y_a = R * np.sin(angle_a)
	# cartesian coordinates of endpoint b
	x_b = R * np.cos(angle_b); y_b = R * np.sin(angle_b)

	chord_1 = patches.PathPatch(patches.Path([[x_a, y_a], [x_b, y_b]],
        [patches.Path.MOVETO, patches.Path.LINETO]),
		edgecolor = color, lw = width, fill = False)

	return ax.add_patch(chord_1)

def get_triangle_side_length(R):
    return R * np.linalg.norm([np.cos(2* np.pi / 3) - 1, np.sin(2* np.pi / 3)])

def l2dist(array_a, array_b):
    assert(len(array_a) == len(array_b))
    
    diffs = [np.linalg.norm([array_a[i][0] - array_b[i][0],
                            array_a[i][1] - array_b[i][1]]) \
                            for i in range(len(array_a))]
    return np.array(diffs)


####

def plot_chord_length_distribution(ax, lengths, R, solution_pct,
												triangle_edge_length):

    # plot distribution of chord lengths
    _n, _bins, _patches = ax.hist(lengths, bins = 100, range = (0, 2*R),
                                     cumulative = True, normed = True)
    fraction_below = round(1 - solution_pct/100, 3)
    _tmp = ax.set(title = "Chord Length Distribution",
    	xlabel = 'Chord Length', ylabel = 'Cumulative Fraction',
    	ylim = (0, max(_n)), yticks = [0.0, fraction_below, 1.0])
    _tmp = ax.plot((triangle_edge_length, triangle_edge_length), (0, max(_n)),
                      c = "red", linewidth = 3)
    _tmp = ax.plot((0, 2*R), (fraction_below, fraction_below), 
                                      color = 'black', linewidth = 2)
    return _tmp

def plot_sample_chords(ax, solution_coord_a, solution_coord_b, R):

	M = len(solution_coord_b)
	assert(M == len(solution_coord_b))

	# construct the chords from the coordinates of their endpoints
	chords_A = [patches.PathPatch(patches.Path([solution_coord_a[i],
		solution_coord_b[i]], [patches.Path.MOVETO, patches.Path.LINETO]),
		edgecolor = "blue", lw = 0.5, fill = False) for i in range(M)]

	_tmp = ax.set(xlim = (-R - 0.01, R + 0.01),
		ylim = (-R - 0.01, R + 0.01), aspect = 1)

	# plot the chords
	for chord in chords_A:
	    _tmp = ax.add_patch(chord)
	return _tmp

def plot_triangle(ax, R):
	triangle = patches.RegularPolygon((0,0), 3, R, lw = 3, fill = False)
	_tmp = ax.add_patch(triangle)
	return _tmp

def plot_circle(ax, c, r, alpha):
	circle = patches.Circle(c, r, ec = 'black', fc = "grey", alpha = alpha)
	_tmp = ax.add_patch(circle)
	return _tmp

def plot_length_distr_and_chords(lengths, R, solution_pct,
			triangle_edge_length, solution_coord_a, solution_coord_b, N):
	fig, ax = plt.subplots(1, 2, figsize = (12, 5))
	ax[1].axis('off')

	# plot distribution of chord lengths
	_tmp = plot_chord_length_distribution(ax[0], lengths, R,
							solution_pct, triangle_edge_length)


	# plot sample of M (out of N) chords
	M = min(500, N)
	_tmp = plot_sample_chords(ax[1], solution_coord_a[:M],
									solution_coord_b[:M], R)
	_tmp = plot_triangle(ax[1], R)

	return

def plot_chords_and_small_circle(c, r, R, N, sol_A_coord_a, sol_A_coord_b,
								sol_B_coord_a, sol_B_coord_b,
								sol_C_coord_a, sol_C_coord_b):
	fig, ax = plt.subplots(1, 3, figsize = (15, 5))
	for i in range(len(ax)): ax[i].axis('off')

	M = min(500, N)
	_tmp = plot_sample_chords(ax[0], sol_A_coord_a[:M], sol_A_coord_b[:M], R)
	_tmp = plot_circle(ax[0], c, r, 0.75)

	_tmp = plot_sample_chords(ax[1], sol_B_coord_a[:M], sol_B_coord_b[:M], R)
	_tmp = plot_circle(ax[1], c, r, 0.75)

	_tmp = plot_sample_chords(ax[2], sol_C_coord_a[:M], sol_C_coord_b[:M], R)
	_tmp = plot_circle(ax[2], c, r, 0.75)

	return

def plot_solution_for_small_vs_large(solution_A_pct, sol_A_small_pct,
	solution_B_pct, sol_B_small_pct, solution_C_pct, sol_C_small_pct):

	fig, ax = plt.subplots(1, 3, figsize = (15, 5))

	_tmp = ax[0].set(xlim = (-0.5, 1.5), ylim = (0,100), xticks = [0, 1],
	                 xticklabels = ['Large', 'Small'],
	                 yticks = [0, int(solution_A_pct), 100])
	_tmp = ax[0].bar([0, 1], [solution_A_pct, sol_A_small_pct], width = 0.5,
	                align = 'center', color = "blue")

	_tmp = ax[1].set(xlim = (-0.5, 1.5), ylim = (0,100), xticks = [0, 1],
	                 xticklabels = ['Large', 'Small'], 
	                 yticks = [0, int(solution_B_pct), 100])
	_tmp = ax[1].bar([0, 1], [solution_B_pct, sol_B_small_pct], width = 0.5,
	                align = 'center', color = "red")
	                 
	_tmp = ax[2].set(xlim = (-0.5, 1.5), ylim = (0,100), xticks = [0, 1],
	                 xticklabels = ['Large', 'Small'], 
	                 yticks = [0, int(solution_C_pct), 100])
	_tmp = ax[2].bar([0, 1], [solution_C_pct, sol_C_small_pct], width = 0.5,
	                align = 'center', color = "green")

	return

####

def point_to_line_dist(coord_p, coord_a, coord_b):
    vector_1 = coord_a - coord_b
    vector_2 = coord_p - coord_b
    len_1 = np.linalg.norm(vector_1)
    len_2 = np.linalg.norm(vector_2)
    angle = np.pi/2 - np.arccos(np.sum(vector_1 * vector_2) / (len_1 * len_2))
    dist = len_2 * np.cos(angle)
    return dist

def get_chord_length(midpoint_radius, circle_radius):
    R_sq = np.square(circle_radius)
    r_sq = np.square(midpoint_radius)
    return 2 * np.sqrt(R_sq - r_sq)

def get_small_lengths(c, r, coord_a_list, coord_b_list):
    small_lengths = []
    for i in range(len(coord_a_list)):
        coord_a = np.array(coord_a_list[i])
        coord_b = np.array(coord_b_list[i])
        dist = point_to_line_dist(c, coord_a, coord_b)
        if dist <= r:
            chord_length = get_chord_length(dist, r)
            small_lengths.append(chord_length)
    return np.array(small_lengths)