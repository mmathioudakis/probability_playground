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

def plot_sample_chords(ax, solution_coord_a, solution_coord_b):

	M = len(solution_coord_b)
	assert(M == len(solution_coord_b))
	# construct the chords from the coordinates of their endpoints
	chords_A = [patches.PathPatch(patches.Path([solution_coord_a[i],
		solution_coord_b[i]], [patches.Path.MOVETO, patches.Path.LINETO]),
		edgecolor = "blue", lw = 0.5, fill = False) for i in range(M)]

	_tmp = ax.set(xlim = (-1.1, 1.1), ylim = (-1.1, 1.1), aspect = 1)
	_tmp = ax.axis('off')

	# plot the triangle and the chords
	for chord in chords_A:
	    _tmp = ax.add_patch(chord)
	triangle = patches.RegularPolygon((0,0), 3, 1, lw = 3, fill = False)
	_tmp = ax.add_patch(triangle)
	return _tmp

def plot_length_distr_and_chords(lengths, R, solution_pct,
			triangle_edge_length, solution_coord_a, solution_coord_b, N):
	fig, ax = plt.subplots(1, 2, figsize = (12, 5))

	# plot distribution of chord lengths
	_tmp = plot_chord_length_distribution(ax[0], lengths, R,
							solution_pct, triangle_edge_length)


	# plot sample of M (out of N) chords
	M = min(500, N)
	_tmp = plot_sample_chords(ax[1], solution_coord_a[:M],
										solution_coord_b[:M])

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