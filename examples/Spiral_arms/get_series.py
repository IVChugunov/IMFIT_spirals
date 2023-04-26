import argparse
from scipy.optimize import minimize, differential_evolution
import numpy as np
from matplotlib import pyplot as plt
import warnings

warnings.filterwarnings("ignore")
default_color = 'green'


class InputParameters(object):
    def __init__(self, is_clockwise, cutoff, points, color, x0=None, y0=None, inc=None, pa=None):
        self.is_clockwise = is_clockwise
        self.cutoff = cutoff
        self.points = points
        self.x0 = x0
        self.y0 = y0
        self.inc = inc
        self.pa = pa
        self.color = color


class SpiralParameters(object):
    def __init__(self, m0, m1, m2, m3, pa, inc, fi0, r0, fi_max, n, color):
        self.m0 = m0
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.pa = pa
        self.inc = inc
        self.fi0 = fi0
        self.r0 = r0
        self.fi_max = fi_max
        self.n = n
        self.color = color


class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __repr__(self):
        return 'x : ' + str(self.x) + ' y : ' + str(self.y)

    def __str__(self):
        return 'x : ' + str(self.x) + ' y : ' + str(self.y)


class PolarPoint(object):
    def __init__(self, fi, r):
        self.fi = fi
        self.r = r

    def __repr__(self):
        return 'fi : ' + str(self.fi) + ' r : ' + str(self.r)


def read_data(filename, parameters_parsers):
    data = open(filename).readlines()
    par = dict()
    spirals = dict()

    names = list(parameters_parsers.keys())
    for line in data:
        parts = line.split()
        if parts[0] in names:
            if parts[0] == 'cutoff':
                if len(parts) == 3:
                    spirals[parts[1]] = parameters_parsers[parts[0]](parts[2])
                else:
                    spirals[default_color] = parameters_parsers[parts[0]](parts[1])
            else:
                par[parts[0]] = parameters_parsers[parts[0]](parts[1])
        else:
            raise NameError(fr'unknown argument name: {parts[0]}')
    return par, spirals


def read_regions(filename):
    data = open(filename).readlines()
    sp = dict()
    for line in data:
        parts = line.split()
        color = default_color
        point = None
        if 'color' in parts[-1]:
            color = parts[-1].split('=')[1]

        if 'point(' in line:
            coords = parts[0].split(',')
            point = Point(float(coords[0][6:]), float(coords[1][:-1]))
        if 'circle(' in line:
            coords = parts[0].split(',')
            point = Point(float(coords[0][7:]), float(coords[1]))

        if point is None:
            continue
        if color not in sp.keys():
            sp[color] = []

        sp[color].append(point)
    return sp


def do_transform(p, pa, inc):
    return [Point(i.x * np.cos(pa) + i.y * np.sin(pa), (- i.x * np.sin(pa) + i.y * np.cos(pa)) / np.cos(inc))
            for i in p]


def to_polar(p, is_clockwise):
    fi00 = np.arctan2(p[0].y, p[0].x)
    if fi00 < 0:
        fi00 += np.pi * 2

    if is_clockwise:
        f_fi0 = np.pi - fi00
    else:
        if fi00 <= np.pi:
            f_fi0 = fi00
        else:
            f_fi0 = fi00 - 2 * np.pi

    r = [np.sqrt(p[i].x ** 2 + p[i].y ** 2) for i in range(len(p))]

    fi = [np.arctan2(p[i].y, -p[i].x if is_clockwise else p[i].x) - f_fi0 for i in range(len(p))]
    fi[0] = 0
    fi = [i if i >= 0 else i + np.pi * 2 for i in fi]
    for i in range(len(fi) - 1):
        if fi[i + 1] < fi[i]:
            fi[i + 1:] = [k + 2 * np.pi for k in fi[i + 1:]]
    return [PolarPoint(fi[i], r[i]) for i in range(len(p))]


def get_function_for_optimization():
    def function_for_optimization(par, *data):
        j = 0

        x0, j = if_null(data[0].x0, par[j], j)
        y0, j = if_null(data[0].y0, par[j], j)

        pa, j = if_null(data[0].pa, par[j], j)
        inc, j = if_null(data[0].inc, par[j], j)

        zero_point = Point(x0, y0)

        div = 0
        for i in range(len(data)):
            points = data[i].points
            points_count = len(points)
            centring_points = [i - zero_point for i in points]
            t_points = do_transform(centring_points, pa, inc)
            polar = to_polar(t_points, data[i].is_clockwise)
            ufi = np.array([k.fi / (2 * np.pi) for k in polar])
            tg_p_angles = par[j + i * 4] + par[j + 1 + i * 4] * ufi + par[j + 2 + i * 4] * ufi ** 2 + par[j + 3 + i * 4] * ufi ** 3
            exps = np.exp([tg_p_angles[i] * polar[i].fi for i in range(len(polar))])
            r = [(polar[0].r * exps[i] - polar[i].r) ** 2 * (1 - args.rate * i / points_count) for i in range(len(exps))]
            div += sum(r)
        return div

    return function_for_optimization


def get_approx_parameters(data, de, acc):
    approx_function = get_function_for_optimization()

    bounds = []

    # bounds for center coordinates
    if data[0].x0 is None:
        bounds.append((min([i.x for i in data[0].points]), max([i.x for i in data[0].points])))

    if data[0].y0 is None:
        bounds.append((min([i.y for i in data[0].points]), max([i.y for i in data[0].points])))

    # pa
    if data[0].pa is None:
        bounds.append((0, 2 * np.pi))

    # inc
    if data[0].inc is None:
        bounds.append((0, np.pi / 2))

    for i in range(len(data)):
        bounds.extend([(0, np.pi / 2), (-1, 1), (-1, 1), (-1, 1)])

    if de:
        last_res = differential_evolution(approx_function, bounds, args=tuple(data)).x
        print("LSA : " + str(approx_function(last_res, *data)))
    else:
        last_res = [(i[0] + i[1]) / 2 for i in bounds]

    res = minimize(approx_function, last_res, args=tuple(data), method='Nelder-Mead', bounds=bounds).x
    print("LSA : " + str(approx_function(res, *data)))

    methods = ['Nelder-Mead']
    while abs(approx_function(last_res, *data) - approx_function(res, *data)) > acc:
        last_res = res
        for method in methods:
            res = minimize(approx_function, res, args=tuple(data), method=method, bounds=bounds).x
            print("LSA : " + str(approx_function(res, *data)))

    j = 0
    x0, j = if_null(data[0].x0, res[j], j)
    y0, j = if_null(data[0].y0, res[j], j)

    zero_point = Point(x0, y0)

    pa, j = if_null(data[0].pa, res[j], j)
    inc, j = if_null(data[0].inc, res[j], j)

    spirals = []
    for i in range(len(data)):
        points = [i - zero_point for i in data[i].points]
        t_points = do_transform(points, pa, inc)

        fi0 = np.arctan2(t_points[0].y, t_points[0].x)

        polar_points = to_polar(t_points, data[i].is_clockwise)

        fi0 = fi0 if fi0 > 0 else fi0 + np.pi * 2
        r0 = polar_points[0].r
        fi_max = polar_points[-1 - data[i].cutoff].fi
        n = polar_points[-1].fi - polar_points[-1 - data[i].cutoff].fi
        spirals.append((polar_points, SpiralParameters(res[j + i * 4], res[j + 1 + i * 4],
                                                       res[j + 2 + i * 4], res[j + 3 + i * 4],
                                                       pa, inc, fi0, r0, fi_max, n, data[i].color)))
    return spirals, zero_point


def if_null(nullable, val, ind):
    if nullable is None:
        return val, ind + 1
    return nullable, ind


def main():
    pars = {'is_clockwise': lambda a: a == "true", 'cutoff': lambda a: int(a), 'X0': lambda a: float(a),
            'Y0': lambda a: float(a), 'inc': lambda a: float(a) * np.pi / 180, 'PA': lambda a: float(a) * np.pi / 180}

    input_data = []
    settings, spirals = read_data(args.file, pars)
    input_points = read_regions(args.regfile)
    for spiral in spirals.keys():
        input_pars = InputParameters(settings['is_clockwise'], spirals[spiral], input_points[spiral], spiral)
        if 'X0' in settings.keys():
            input_pars.x0 = settings['X0']

        if 'Y0' in settings.keys():
            input_pars.y0 = settings['Y0']

        if 'inc' in settings.keys():
            input_pars.inc = settings['inc']

        if 'PA' in settings.keys():
            input_pars.pa = settings['PA']

        input_data.append(input_pars)

    result, p0 = get_approx_parameters(input_data, args.de, args.acc)

    for j in range(len(result)):
        theta = np.array([i for i in np.arange(0, result[j][0][-1].fi + np.pi / 8, 0.01)])
        ufi = theta / (2 * np.pi)
        tg_pitch_angles = result[j][1].m0 + result[j][1].m1 * ufi + result[j][1].m2 * ufi ** 2 + result[j][
            1].m3 * ufi ** 3
        r_model = result[j][1].r0 * np.exp([theta[i] * tg_pitch_angles[i] for i in range(len(theta))])
        pitch_angles = np.arctan(
            result[j][1].m0 + 2 * result[j][1].m1 * ufi + 3 * result[j][1].m2 * ufi ** 2 + 4 * result[j][
                1].m3 * ufi ** 3)

        ax1 = plt.subplot(1, 2, 1, projection='polar')
        ax1.set_rlim(0, max([i.r for i in result[j][0]]) * 1.2)
        ax1.plot([i.fi for i in result[j][0]], [i.r for i in result[j][0]], 'rx' if args.x else 'y-')
        ax1.plot(theta, r_model)

        ax2 = plt.subplot(1, 2, 2)
        ax2.set_xlim(0, max([i.r for i in result[j][0]]) * 1.2)
        pitch_angles = [i * 180 / np.pi for i in pitch_angles]
        ax2.plot(r_model[0:len(pitch_angles)], pitch_angles[0:])
        ax2.set_yticks(np.arange(min(pitch_angles), max(pitch_angles), (max(pitch_angles) - min(pitch_angles)) / 10))
        plt.savefig(result[j][1].color + '.pdf')

        plt.cla()
        plt.clf()

        lines = [
            'x0 ' + str(p0.x) + ' fixed\n',
            'y0 ' + str(p0.y) + ' fixed\n',
            'm0 ' + str(result[j][1].m0) + ' fixed\n',
            'm1 ' + str(result[j][1].m1) + ' fixed\n',
            'm2 ' + str(result[j][1].m2) + ' fixed\n',
            'm3 ' + str(result[j][1].m3) + ' fixed\n',
            'PA ' + str(result[j][1].pa * 180 / np.pi) + ' fixed\n',
            'inc ' + str(result[j][1].inc * 180 / np.pi) + ' fixed\n',
            'fi0 ' + str(result[j][1].fi0 * 180 / np.pi) + ' fixed\n',
            'r0  ' + str(result[j][1].r0) + ' fixed\n',
            'fi_max ' + str(result[j][1].fi_max * 180 / np.pi) + ' fixed\n',
            'n ' + str(result[j][1].n * 180 / np.pi) + ' fixed\n'
        ]
        open(result[j][1].color + '_parameters.dat', 'w').writelines(lines)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate series coefficients')
    parser.add_argument('file', metavar='filename', type=str, help='file with setting')
    parser.add_argument('regfile', metavar='regfilename', type=str, help='file with points in reg format')
    parser.add_argument('-de', action="store_true", help='Add differential evolution method')
    parser.add_argument('-x', action="store_true", help='Add crosses to plot')
    parser.add_argument('-acc', type=float, default=0.1, help='Accuracy')
    parser.add_argument('-rate', type=float, default=0, help='Weights decrease rate')

    args = parser.parse_args()
    main()
