#!/usr/bin/env python
# coding: utf-8

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.decomposition import PCA
import powerlaw
from scipy.optimize import curve_fit
from scipy.stats import powerlaw as pw
from scipy.stats import lognorm, expon
from scipy.stats import kstest, ks_2samp, anderson_ksamp
import math
import csv
import cv2
import shapefile
import sys
from time import time
from multiprocessing import Process, Array, Lock


THETA = math.pi/180

def r2_score(y_true, y_predicted):
    sse = sum((y_true - y_predicted)**2)
    tse = (len(y_true) - 1) * np.var(y_true, ddof=1)
    r2 = 1 - (sse / tse)
    return r2


###############################################################################
# File handling
###############################################################################


def load_lines_from_file(File, offsetX, offsetY, k=0, segm_groups=[]):

    # Read the file header and identify the drawn lines start and finish
    with open(File) as csv_file:
       csv_reader = csv.reader(csv_file, delimiter=' ')
       row_count = 0
       for row in csv_reader:
           for cell in row:
               if cell == 'FreeHand':
                   lines_start = int(row[1]) - 1
                   lines_end = int(row[2]) + int(row[1]) + 1
                   break

    lines = []
    segment = []
    with open(File) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        row_count = 0
        segment = []
        for row in csv_reader:
            segment = []
            if row_count < lines_end and row_count > lines_start:
                collumn_count = 0
                for collumn in row:
                    # Ignore collumn indexes before line 3
                    if collumn_count > 2:
                        try:
                            # Get the value in each collumn and replace the
                            # commas for dots
                            value = float(collumn.replace(",", "."))
                            segment.append(value)
                        except ValueError:
                            print("", end="")
                    collumn_count = collumn_count + 1
                try:
                    segment = np.reshape(segment, (int(np.size(segment)/3), 3))
                except Exception as e:
                    print(e)
                    break
            
            segm_group_row = []

            # Get individual segments
            for i in range(0, np.shape(segment)[0]-1):
                lines.append([segment[i][0], -segment[i][2], segment[i+1][0],
                            -segment[i+1][2]])
                segm_group_row.append(k)
                k = k + 1

            if np.size(segm_group_row) > 0:
                segm_groups.append(segm_group_row)
            row_count = row_count + 1
            
    lines = np.reshape(lines, (np.shape(lines)[0], 4))

    minX = min(np.append(lines[:, 0], lines[:, 2]))
    minY = min(np.append(lines[:, 1], lines[:, 3]))
    maxX = max(np.append(lines[:, 0], lines[:, 2]))
    maxY = max(np.append(lines[:, 1], lines[:, 3]))

    for i in range(0, np.shape(lines)[0]):
        lines[i][0] = (lines[i][0]-minX+offsetX)
        lines[i][1] = (lines[i][1]-minY+offsetY)
        lines[i][2] = (lines[i][2]-minX+offsetX)
        lines[i][3] = (lines[i][3]-minY+offsetY)

    maxY = max(np.append(lines[:, 1], lines[:, 3]))
    maxX = max(np.append(lines[:, 0], lines[:, 2]))

    return lines, maxX, maxY, segm_groups


def load_lines_from_xpp_file(File, offsetX=0, offsetY=0, k=0, segm_groups=[]):
  #File = "Gaivota2D.xpp"
  lines = []
  segm_groups = []
  points = []
  x0 = 0
  y0 = 0
  with open(File) as csv_file:
      points = []
      csv_reader = csv.reader(csv_file, delimiter='|')

      row_count = 0
      for row in csv_reader:
        for cell in row:
          if cell == 'FreeHand':
            planes_start = int(row[2])
            planes_end = int(row[2]) + int(row[1])
            break
        if 'planes_start' in locals():
            if row_count >= planes_start and row_count < planes_end:
              points = []

              row_values = row[6].split()


              for i in range(np.size(row_values)):
                if i > 2 and i < np.size(row_values):
                  points.append(float(row_values[i]))
              points = np.reshape(points, (int(np.size(points)/3), 3))
              points = np.asarray([points[:,0], points[:,2], -points[:,1]]).T
              #print(points)
        
              X = points[:,0]
              Y = -points[:,1]

              pca = PCA(n_components=1)
              data = np.transpose(np.asarray((X, Y)))
              pca.fit(data)
              X_pca = pca.transform(data)
              X_new = pca.inverse_transform(X_pca)

              # Obter os extremos horizontais
              p = np.where(X_new[:, 0] == np.amax(X_new[:, 0]))[0][0]
              q = np.where(X_new[:, 0] == np.amin(X_new[:, 0]))[0][0]

              # Verificar se a linha é vertical para calcular os extremos verticais
              if p == q:
                  p = np.where(X_new[:, 1] == np.amax(X_new[:, 1]))[0][0]
                  q = np.where(X_new[:, 1] == np.amin(X_new[:, 1]))[0][0]

              line = [X_new[p][0] + x0, X_new[p][1] + y0,
                      X_new[q][0] + x0, X_new[q][1] + y0]

              lines.append(line)
              #print(line)
        row_count +=1

  return lines


def save_shapefile(file_address, lines, angles):
    w = shapefile.Writer(file_address, shapeType=3)
    w.field('fracture_i', 'N')
    w.field('fractdir', 'N')
    w.field('fractlength', 'N')

    for i in range(0, np.shape(lines)[0]):
        point0 = lines[i][0], lines[i][1]
        point1 = lines[i][2], lines[i][3]
        # Add record
        w.record(i, angles[i, 0], angles[i, 1])
        # Add geometry
        w.line([[[point0[0], point0[1]], [point1[0], point1[1]]]])

    w.close()
    return


###############################################################################
# Geometry
###############################################################################


def point_distance(a, b):
    x0, y0 = a
    x1, y1 = b
    d = math.sqrt(math.pow((x0 - x1), 2) + math.pow((y0 - y1), 2))
    return d

def compute_distance(x0, y0, x1, y1):
    #dist = math.sqrt(math.pow((x1 - x0), 2) + math.pow((y1 - y0), 2))
    dist = (math.pow((x1 - x0), 2) + math.pow((y1 - y0), 2))**.5
    return dist

def get_intersection_coordinates(circle_center, radius, line_start, line_end):
    dx = line_end[0] - line_start[0]
    dy = line_end[1] - line_start[1]
    a = dx**2 + dy**2
    b = 2 * (dx * (line_start[0] - circle_center[0]) + dy * (line_start[1] - circle_center[1]))
    c = (line_start[0] - circle_center[0])**2 + (line_start[1] - circle_center[1])**2 - radius**2
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        return []
    elif discriminant == 0:
        u = -b / (2 * a)
        return [(line_start[0] + u * dx, line_start[1] + u * dy)]
    else:
        u1 = (-b + discriminant**0.5) / (2 * a)
        u2 = (-b - discriminant**0.5) / (2 * a)
        return [(line_start[0] + u1 * dx, line_start[1] + u1 * dy), (line_start[0] + u2 * dx, line_start[1] + u2 * dy)]


def circular_scanline(lines, circle_center, radius):
    segments_in_circle = []
    for i in range(np.shape(lines)[0]): # np.shape(lines)[0]
        line = lines[i]
        intersection = get_intersection_coordinates(circle_center, radius, (line[0], line[1]), (line[2], line[3]))
        if np.shape(intersection)[0] == 2:
            a, b = intersection
            x1, y1 = a
            x2, y2 = b

            dist_a = point_distance(circle_center, (line[0], line[1]))
            dist_b = point_distance(circle_center, (line[2], line[3]))
            
            if dist_a <= radius and dist_b <= radius:
                segments_in_circle.append([line[0], line[1], line[2], line[3]])
                
            if dist_a > radius and dist_b > radius:
                dist1 = point_distance((line[0], line[1]), (x1, y1)) + point_distance((x1, y1), (line[2], line[3]))
                dist2 = point_distance((line[0], line[1]), (line[2], line[3]))
                if math.isclose(dist1, dist2):
                    segments_in_circle.append([x1, y1, x2, y2])
                   
            if dist_a <= radius and dist_b > radius:
                if point_distance(a, (line[2], line[3])) < point_distance(b, (line[2], line[3])):
                    segments_in_circle.append([line[0], line[1], x1, y1])
                else:
                    segments_in_circle.append([line[0], line[1], x2, y2])
                    
            if dist_a > radius and dist_b <= radius:
                if point_distance(a, (line[0], line[1])) < point_distance(b, (line[0], line[1])):
                    segments_in_circle.append([line[2], line[3], x1, y1])
                else:
                    segments_in_circle.append([line[2], line[3], x2, y2])

    segments_in_circle = np.asarray(segments_in_circle)
    # print(np.shape(segments_in_circle))
    return segments_in_circle

def compute_line_angles(lines):
    # Get distance and angle of all lines regarding north
    n = np.shape(lines)[0]
    angles = np.zeros((n, 2), np.float64)
    for i in range(0, n):
        length = compute_distance(lines[i][0], lines[i][1], lines[i][2],
                                  lines[i][3])

        if lines[i][1] == lines[i][3]:
            angle = 90
        else:
            angle = (math.atan((lines[i][2]-lines[i][0])
                               / (lines[i][1]-lines[i][3])) / (THETA))

        if (angle < 0):
            angle = angle + 360
        if (angle > 180):
            angle = angle - 180
        if math.isnan(angle):
            angle = 0

        angles[i, 0] = np.uint(angle)
        angles[i, 1] = length
    return angles


def transfom_to_geocoordinates(lines, base_offset):
    lines[:,0] = base_offset[0]+lines[:,0]
    lines[:,1] = base_offset[1]-lines[:,1]
    lines[:,2] = base_offset[0]+lines[:,2]
    lines[:,3] = base_offset[1]-lines[:,3]
    return lines

def transfom_to_geocoordinates_point(point, base_offset):
    point[0] = base_offset[0]+point[0]
    point[1] = base_offset[1]-point[1]
    return point

###############################################################################
# Functions to estimate power law distributions 
###############################################################################

def powerlaw_MLE(data): # Rizzo, 2017. pag 19 Clauset, 2009
    # Standard MLE for powerlaw distributions
    # TODO: Use SGD to get optmize value for fractal dimension alpha (D) 
    alpha = 1 + np.size(data)*np.power(abs(np.sum(np.log(data/np.min(data)))), -1)
    return alpha


def powerlaw_fit_MLE(data):
  # Method proposed by Deluca et al, 2013 to estimate alpha and xmin attributes in
  # powerlaw distributions.
  # data = np.copy(set2_4[:,1])
  distances = []
  alphas = []

  data = np.sort(data)

  for i in range(int(np.size(data)*0.75)):

    temp_data = data[i:]
    x = data[i]
    if np.size(temp_data) > 2:

      alpha_mle = powerlaw_MLE(temp_data)
      alphas.append(alpha_mle)

      # fit = powerlaw.Fit(temp_data, x)
      # alpha_mle = fit.alpha
      # alphas.append(alpha_mle)

      x_ = np.linspace(pw.ppf(0.01, alpha_mle),
                    pw.ppf(0.99, alpha_mle), np.size(temp_data))

      P = x*np.power((1-x_),(-1/(alpha_mle-1)))
      D = ks_2samp(temp_data, P)[0]

      ks_t = kstest_hit_ratio2(P, alpha_mle, x, 1000, D, 'distance') # S or D
    #   ks_t = kstest_hit_ratio2(temp_data, alpha_mle, x, 1000)
      distances.append(ks_t)

  xmin = data[np.where(distances == np.max(distances))][0]
  alpha = np.asarray(alphas)[np.where(distances == np.max(distances))][0]
  return alpha, xmin


###############################################################################
# Plot functions 
###############################################################################

def rosechart(angles_, weights=None, filename=None):
  bin_edges = np.arange(-5, 366, 10)
  number_of_strikes, bin_edges = np.histogram(angles_, bin_edges, weights)

  number_of_strikes[0] += number_of_strikes[-1]
  half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
  two_halves = np.concatenate([half, half])

  fig = plt.figure() # figsize=(10, 10)
  ax = fig.add_subplot(111, projection='polar')
  ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves,
          width=np.deg2rad(10), bottom=0.0, color='gray', edgecolor='lightgray')
  ax.set_theta_zero_location('N')
  ax.set_theta_direction(-1)
  ax.set_thetagrids(np.arange(0, 360, 90), labels=['','','',''])
  #ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
  # ax.set_rgrids(np.arange(1, two_halves.max()+1, 2),angle=0,weight='black')
  ax.set_title("n = "+str(np.size(angles_)), y=1, fontsize=14)
  ax.set_yticks([])
  fig.canvas.draw()

  # Save it to a numpy array.
  plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
  plot = plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

  if filename is not None:
    plt.savefig(filename)

  #plt.show()
  return plot


def draw_lines(canvas, lines, x_offset, y_offset, color):
    for i in range(np.shape(lines)[0]): #np.shape(lines)[0]
        line = lines[i]
        canvas = cv2.line(canvas, (int(line[0]+ x_offset), int(line[1]+ y_offset)),
                          (int(line[2]+ x_offset), int(line[3]+ y_offset)), color)
    return canvas


def plot_history(n_seeds, cost_history, max_iterations, filename=None):
  fig, ax = plt.subplots(figsize=(6, 4))
  for i in range(n_seeds):
    ax.plot(cost_history[i])
  ax.set_xlabel('Iterations', fontsize=14)
  ax.set_ylabel('Cost', fontsize=14)

  if filename is not None:
    plt.savefig(filename)

  plt.show()


def plot_powerlaw_fitting(X=[], alpha_=1.1, xmin_=1, plot=True, filename=None):
  X = X[X > xmin_]
  X = X.flatten()
  X = np.sort(X)
  X = np.flip(X)

  S = 1-powerlaw_gen_cdf(alpha_, np.min(X), np.max(X), np.size(X))
  A2 = anderson_ksamp([get_cdf(X), S], midrank=True)[0]
  HD = hellinger(get_cdf(X), S)
  ks_d = ks_2samp(get_cdf(X), S)[0]
  ks_t = kstest_hit_ratio2(X, alpha_, xmin_)

  # x_ = np.linspace(pw.ppf(0.01, alpha_),
  #                   pw.ppf(0.99, alpha_), np.size(X))
  # P = xmin_*np.power((1-x_),(-1/(alpha_-1)))
  # # S = temp_data/np.max(temp_data)
  # D = ks_2samp(X, P)[0]
  # ks_t = kstest_hit_ratio2(X, alpha_, xmin_, ref=D, mode='distance')

  uniform_data = np.linspace(xmin_, np.max(X), np.size(X), dtype=np.float64)
  Y = 1-(uniform_data/xmin_)**(-alpha_+1)  # Equation 2.6 from Clauset, 2009
  # Y = generate_powerlaw_data(alpha, xmin_, np.size(X), random=False)

  Y_ = np.copy(Y)
  Y = Y.cumsum()
  Y = Y/np.max(Y)

  pars, cov = curve_fit(f=power_law, xdata=X, ydata=Y, p0=[0, 0], bounds=(-np.inf, np.inf))
  r2 = r2_score(Y, power_law(X, *pars))


  if plot==True:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.grid(True)
    # ax.plot(X, Y, 'bx', label="Data")
    ax.plot(X, Y, '.', color='#00b3b3', label='Data')

    ax.plot(X, power_law(X, *pars), linewidth=2, linestyle='--', color='black', label='Power-law fitting')
    ax.plot([], [], ' ', label=r'$\hat\alpha$: ' + str("%.1f" % alpha_))
    ax.plot([], [], ' ', label=r'$x_{min}$: ' + str("%.1f" % xmin_))
    
    ax.legend(loc='lower left',fontsize=16)
    # ax.set_title(r'$\alpha$: ' + str("%.1f" % alpha) + r'     $x_{min}$: ' + str("%.1f" % xmin), fontsize=14)
    ax.set_xlabel(r'$x$', fontsize=18)
    ax.set_ylabel(r'$P(x)$', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    ax.set_yscale('log')
    ax.set_xscale('log')

    # Edit the major and minor ticks of the x and y axes
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')

    # Edit the major and minor tick locations of x and y axes
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))

    if filename is not None:
      plt.savefig(filename)

    plt.show()

  return r2, ks_t, ks_d, A2, HD
# ref: https://towardsdatascience.com/basic-curve-fitting-of-scientific-data-with-python-9592244a2509


def plot_distributions_fitting(data, xmin, filename=None, plot=True):
    # def get_cdf(data):
    #     # return np.histogram(data, np.size(data), density=False)[0]
    #     H,X1 = np.histogram(data, np.size(data), density = True )
    #     dx = X1[1] - X1[0]
    #     F1 = np.cumsum(H)*dx
    #     return F1
    
    # def powerlaw_gen_cdf(alpha, xmin, xmax, size):
    #     uniform_data = np.linspace(xmin, xmax, size, dtype=np.float64)
    #     # Equation 2.6 in Clauset et al. 2009
    #     Y = 1-(uniform_data/xmin)**(-alpha+1)
    #     return Y

    X = np.copy(data)
    xmin_ = xmin
    plot = True

    X = X[X >= xmin_]
    X = X.flatten()
    X = np.sort(X)
    X = np.flip(X)

    original = np.copy(data)
    original = np.flip(np.sort(original.flatten()))
    original = original[original < xmin_]

    fit = powerlaw.Fit(X, xmin=xmin_) # powerlaw
    shape, loc, scale = lognorm.fit(X, floc=0) # lognormal
    fittedlog = lognorm(shape, loc, scale)
    loc, scale = expon.fit(X, floc=xmin_) # exponential
    fittedexpon = expon(loc, scale)

    cf = np.linspace(0, 1-1/np.size(X), np.size(X))
    cf2 = np.linspace(1-1/np.size(X), (1-(1/np.size(X)))+np.size(original)*(1/np.size(X)), np.size(original))

    # powerlaw_gen_cdf(fit.alpha, fit.xmin, np.max(X), np.size(X))
    # fit.alpha, fit.xmin
    # plt.plot(np.flip(fit.power_law.ccdf()))
    # plt.plot(fittedexpon.sf(X))
    # plt.show()

    score_powerlaw_ad = anderson_ksamp([cf, np.flip(fit.power_law.ccdf())], midrank=True)[0]
#     if score_powerlaw_ad < 0: score_powerlaw_ad = score_powerlaw_ad + 1
    score_log_ad = anderson_ksamp([cf, fittedlog.sf(X)], midrank=True)[0]
#     if score_log_ad < 0: score_log_ad = score_log_ad + 1
    score_expon_ad = anderson_ksamp([cf, fittedexpon.sf(X)], midrank=True)[0]
#     if score_expon_ad < 0: score_expon_ad = score_expon_ad + 1

    score_powerlaw_hd = hellinger(cf, np.flip(fit.power_law.ccdf()))
    score_log_hd = hellinger(cf, fittedlog.sf(X))
    score_expon_hd = hellinger(cf, fittedexpon.sf(X))

    score_powerlaw_ks = ks_2samp(cf, np.flip(fit.power_law.ccdf()))[0]
    score_log_ks = ks_2samp(cf, fittedlog.sf(X))[0]
    score_expon_ks = ks_2samp(cf, fittedexpon.sf(X))[0]


    fig, ax = plt.subplots(figsize=(6, 6))
    ax.grid(True)
    ax.plot(original, cf2, '.', color='gray')
    ax.axvline(x = xmin_, color = 'gray', linestyle='--', label = 'Cutoff')
    ax.plot(X, cf, '.', color='#7f0000', label='Data')

    fit.power_law.plot_ccdf(ax = ax, color = 'black', linewidth=2, linestyle = '--',
                            label = 'Powerlaw - AD: ' + str(round(score_powerlaw_ad, 3)) +
                                    ', HD: ' + str(round(score_powerlaw_hd, 3)) +
                                    ', KS: ' + str(round(score_powerlaw_ks, 3)))
    # ax.plot(X, power_law(X, *pars), linewidth=2, linestyle='--', color='black', label='Power-law fitting')
    ax.plot(X, fittedlog.sf(X), linewidth=2, linestyle='--', color='blue',
            label='Lognormal - AD: ' + str(round(score_log_ad, 3)) +
                    ', HD: ' + str(round(score_log_hd, 3)) +
                    ', KS: ' + str(round(score_log_ks, 3)) )
    
    ax.plot(X, fittedexpon.sf(X), linewidth=2, linestyle='--', color='red',
            label='Exponential - AD: ' + str(round(score_expon_ad, 3)) + 
                    ', HD: ' + str(round(score_expon_hd, 3)) +
                    ', KS: ' + str(round(score_expon_ks, 3)) )

    # ax.plot([], [], ' ', label=r'$\hat\alpha$: ' + str("%.1f" % alpha))
    # ax.plot([], [], ' ', label=r'$x_{min}$: ' + str("%.1f" % xmin_))
    ax.legend(loc='lower left',fontsize=10)
    # ax.set_title(r'$\alpha$: ' + str("%.1f" % alpha) + r'     $x_{min}$: ' + str("%.1f" % xmin), fontsize=14)
    ax.set_xlabel(r'$x$', fontsize=14)
    ax.set_ylabel(r'$CCDF(x)$', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on') # Edit the major and minor ticks of the x and y axes
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0)) # Edit the major and minor tick locations of x and y axes
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if plot == True:
        plt.show()





###############################################################################
# Metrics of distribution comparison 
###############################################################################


def hellinger (p, q):
  H = (1/np.sqrt(2))*np.sqrt(np.sum(np.power((np.sqrt(p)-np.sqrt(q)), 2)))
  return H

# def ks_distance(X, alpha=1.1, xmin=1):
#   uniform_data = np.linspace(xmin, np.max(X), np.size(X), dtype=np.float64)
#   P = (uniform_data/np.min(uniform_data))**(-alpha+1)

#   X = np.flip(np.sort(X))
#   x_ = X/np.max(X)

#   D = np.max(np.absolute(P - x_))
#   RMSE = np.sqrt(sum(np.power((P - x_), 2)))

#   return D, RMSE


def kstest_hit_ratio2(X, alpha=1.1, xmin=1, max_iter=2500, ref=0.05, mode='pvalue', synthetic_data=None):
  # X = np.sort(X)
  if type(synthetic_data) == type(None):
    # np.random.seed(3)
    # r = np.array(np.random.uniform(0.0, 1.0, size=(max_iter, np.size(X))), dtype=np.float64)
    r = np.random.default_rng().uniform(size=(max_iter, np.size(X)))
    power = lambda r: xmin*np.power((1-r),(-1/(alpha-1)), dtype=np.float64)# Equation D.4 from Clauset, 2009
    g = np.sort(power(r), axis=1)

  temp_X = np.reshape(np.tile(X, max_iter), (max_iter, np.size(X)))
  test = np.asarray(list(map(kstest, temp_X, g)))

  if mode=='pvalue':
    return (np.size(np.where(test[:,1] > ref))/max_iter)*100
  if mode=='distance':
    return (np.size(np.where(test[:,0] > ref))/max_iter)*100


###############################################################################
# Distribution generation
###############################################################################


def powerlaw_gen_data(alpha, xmin, size):
    r = np.random.default_rng().uniform(size=size)
    # Equation D.4 in Clauset et al. 2009
    power = xmin*(1-r)**(-1/(alpha-1))
    return power


def powerlaw_gen_pdf(alpha, xmin, xmax, size):
    uniform_data = np.linspace(xmin, xmax, size, dtype=np.float64)
    # Equation 2.2 in Clauset et al. 2009
    p = ((alpha - 1) / xmin) * (uniform_data / xmin) ** -alpha
    return p


def powerlaw_gen_cdf(alpha, xmin, xmax, size):
    uniform_data = np.linspace(xmin, xmax, size, dtype=np.float64)
    # Equation 2.6 in Clauset et al. 2009
    Y = (uniform_data/xmin)**(-alpha+1)
    return Y


def power_law(x, a, b):
    return a*np.power(x, b)


# def generate_powerlaw_data(alpha, xmin, size, random=False):
#   if random == True:
#     x =np.random.uniform(0.0, 1.0, size=size)
#   else:
#     x = np.linspace(pw.ppf(0.01, alpha),
#                     pw.ppf(0.99, alpha), size)
    
#   P = xmin*np.power((1-x),(-1/(alpha-1)))
#   P = np.sort(P)
#   return P


# def inverse_cumulative_powerlaw(bins, alpha, cumulative=True, x_min=None):
#     # bins = rescale(bins, (x_min, 1))

#     # bins = 1 - bins

#     if x_min == None:
#       x_min = np.min(bins)
#     # power_law = ((alpha-1)/x_min) * (bins/x_min) ** -alpha

#     power_law = (bins/x_min)**(-alpha+1)
#     if cumulative is True:
#         power_law = power_law.cumsum()
#         power_law /= power_law[-1] 
#     return power_law


def get_pdf(data):
  return np.histogram(data, np.size(data), density=True)[0]


def get_cdf(data):
  # return np.histogram(data, np.size(data), density=False)[0]
  H,X1 = np.histogram(data, np.size(data), density = True )
  dx = X1[1] - X1[0]
  F1 = np.cumsum(H)*dx
  return F1


###############################################################################
# SGD 
###############################################################################


def cost_function(data, alpha=None, x_min=None, metric='KS_t'):
  P = np.copy(data[data >= x_min])

  # uniform_data = np.linspace(x_min, np.max(P), np.size(P), dtype=np.float64)
  # S = (uniform_data/np.min(uniform_data))**(-alpha+1)

  # S = 1-powerlaw_gen_cdf(alpha, x_min, np.max(P), np.size(P))
  
  if metric == 'A2':
    S = 1-powerlaw_gen_cdf(alpha, x_min, np.max(P), np.size(P))
    A2 = anderson_ksamp([get_cdf(P), S], midrank=True)[0]
    cost = A2
    
  if metric == 'HD':
    S = 1-powerlaw_gen_cdf(alpha, x_min, np.max(P), np.size(P))
    HD = hellinger(get_cdf(P), S)
    cost = HD

  if metric == 'KS':
    S = 1-powerlaw_gen_cdf(alpha, x_min, np.max(P), np.size(P))
    KS = ks_2samp(get_cdf(P), S)[0]
    print(KS, flush=True)
    cost = KS

  if metric == 'all':
    S = 1-powerlaw_gen_cdf(alpha, x_min, np.max(P), np.size(P))
    A2 = anderson_ksamp([get_cdf(P), S], midrank=True)[0]
    HD = hellinger(get_cdf(P), S)
    KS = ks_2samp(get_cdf(P), S)[0]
    cost = HD + A2 + KS

  if metric == 'KS_t_DH':
    S = 1-powerlaw_gen_cdf(alpha, x_min, np.max(P), np.size(P))
    A2 = anderson_ksamp([get_cdf(P), S], midrank=True)[0]
    HD = hellinger(get_cdf(P), S)
    KS = ks_2samp(get_cdf(P), S)[0]

    ks_t = kstest_hit_ratio2(P, alpha, x_min, max_iter=1000)
    cost = 100-ks_t + HD + A2 + KS

  if metric == 'KS_t':
    ks_t = kstest_hit_ratio2(P, alpha, x_min, max_iter=2500, synthetic_data=None)
    cost = 100-ks_t

  if metric == 'KS_t0':
    x_ = np.linspace(pw.ppf(0.001, alpha),
                    pw.ppf(0.999, alpha), np.size(P))
    S_ = x_min*np.power((1-x_),(-1/(alpha-1)))
    D = ks_2samp(P, S_)[0]
    ks_t = 100-kstest_hit_ratio2(S_, alpha, x_min, 2500, D, 'distance')
    cost = ks_t
    
  return cost


def sgd(data, alpha0, x_min0, learning_rate, momentum, change, fix_x_min, metric):
  alpha1 = alpha0
  x_min1 = x_min0 

  step = ((np.random.rand()*-2)+1)*learning_rate + momentum*change[0]
  alpha1 = alpha0 + step
  change[0] = step
  if alpha1 <= 1.1:
    # alpha1 = alpha0 - step
    # change[0] = -step
    alpha1 = 1.1
    change[0] = 0

  if fix_x_min == None:
    step = (((np.random.rand()*-2)+1)*learning_rate + momentum*change[1])
    x_min1 = x_min0 + step
    change[1] = step
    if x_min1 < np.min(data):
       x_min1 = np.min(data)
       change[1] = 0
    if x_min1 > np.percentile(data, 75):
       x_min1 = np.percentile(data, 75)
       change[1] = 0
    # if x_min1 >= np.min(data) and x_min1 <= np.percentile(data, 75):
    #   x_min1 = x_min0 - step
    #   change[1] = -step
    # else:
    #   x_min1 = x_min0
  else:
    x_min0 = fix_x_min
    x_min1 = fix_x_min

  custo0 = cost_function(data, alpha0, x_min0, metric)
  custo1 = cost_function(data, alpha1, x_min1, metric)

  # alpha_temp = alpha0
  xmin_temp = x_min0

  if custo1 < custo0:
    alpha_temp = alpha1
    alpha0 = alpha1

    # xmin_temp = x_min1
    x_min0 = x_min1
    custo0 = custo1

  return alpha0, x_min0, custo0, change 


# Core function that call the SGD and save the individual results in global arrays
def core_function(data, index, max_iterations, alpha0, x_min0, metric, learning_rate, momentum, fix_x_min, cost_history, alpha_history, x_min_history, early_stopping=True):
  alpha = np.copy(alpha0)
  x_min = np.copy(x_min0)

  index = int(index)

  temp_cost = []
  temp_alpha = []
  temp_xmin = []

  prior_steps = [0,0] # Necessary to compute momentum
  cost_change_history = [] # Necessary for stopping criteria

  count = 0 # Necessary for stopping criteria
  cost0 = 0
  cond = True

  while cond:
    # Sochastic Gradient Descent
    alpha, x_min, cost, prior_steps = sgd(data, alpha, x_min, learning_rate,
                                         momentum, prior_steps, fix_x_min, metric)
    
    temp_cost.append(cost)
    temp_alpha.append(alpha)
    temp_xmin.append(x_min)

    # Stopping criteria
    cost_change_history.append(abs(cost0-cost))
    history_size = np.size(cost_change_history)
    max_stagnation = int(max_iterations*0.1)
    if history_size > max_stagnation and early_stopping:
       if np.size(np.unique(cost_change_history[history_size-max_stagnation:history_size])) == 1:
          cond = False
    if count > max_iterations-2:
       cond = False
    cost0 = cost
    count += 1

  # print(list(cost_change_history.queue))
  cost_history[(index*max_iterations):(index*max_iterations)+np.size(temp_cost)] = temp_cost
  alpha_history[(index*max_iterations):(index*max_iterations)+np.size(temp_cost)] = temp_alpha
  x_min_history[(index*max_iterations):(index*max_iterations)+np.size(temp_cost)] = temp_xmin
    
  return


def powerlaw_fit_SGD(data,max_iterations,learning_rate,n_seeds, momentum=0.8, fix_x_min=None, metric='kS_t', multiprocessing=True, early_stopping=True, history=True):
  # Creating shared variables
  cost_history = Array('d', np.full((n_seeds*max_iterations),-1))
  alpha_history = Array('d', np.full((n_seeds*max_iterations),-1))
  x_min_history = Array('d', np.full((n_seeds*max_iterations),-1))

  # Preparing the input data
  data = np.sort(data)
  data = data.flatten()

  # Generating initial values for alpha and xmin
  alpha0 = 1 + np.random.rand(n_seeds)*2
  x_min0 = np.array(np.random.uniform(low=np.min(data), high=np.percentile(data, 75), size=n_seeds), dtype=np.float64)

  # Run the core fuction that will launch the SGD
  if multiprocessing:
    # start_time = time()
    processes = list()
    for index in range(n_seeds):
      p = Process(target=core_function, args=(data, index, max_iterations, alpha0[index], x_min0[index], metric, learning_rate,
            momentum, fix_x_min, cost_history, alpha_history, x_min_history, early_stopping))
      processes.append(p)
      p.start()
    
    for p in processes:
      p.join()
    # end_time = time()
    # print('Parallel time: '+str(end_time-start_time)+'s')
  else:
    # start_time = time()
    for index in range(n_seeds):
       core_function(data, index, max_iterations, alpha0[index], x_min0[index], metric, learning_rate, momentum,
                     fix_x_min, cost_history, alpha_history, x_min_history, early_stopping)
    # end_time = time()
    # print('Serial execution time: '+str(end_time-start_time)+'s')

  # Preparing the output data
  alpha_history = np.asarray(alpha_history, dtype=np.float64)
  x_min_history = np.asarray(x_min_history, dtype=np.float64)
  cost_history = np.asarray(cost_history, dtype=np.float64)

  alpha_history = np.reshape(alpha_history, (n_seeds, max_iterations))
  x_min_history = np.reshape(x_min_history, (n_seeds, max_iterations))
  cost_history = np.reshape(cost_history, (n_seeds, max_iterations))

  alpha_history[np.where(alpha_history==-1)] = None
  x_min_history[np.where(x_min_history==-1)] = None
  cost_history[np.where(cost_history==-1)] = None

  if history==True:
     return alpha_history, x_min_history, cost_history
  else:
    alpha_sgd = alpha_history[np.where(cost_history==np.nanmin(cost_history))][0]
    x_min_sgd = x_min_history[np.where(cost_history==np.nanmin(cost_history))][0]
    return alpha_sgd, x_min_sgd


if __name__ == '__main__':
    learning_rate = 0.05
    momentum=0.2
 
    
    ## Generating synthetic distribution
    alpha = 2.5
    xmin = 1
    # np.random.seed(0)
    r = np.random.uniform(0.0, 1.0, size=100)
    power = xmin*(1-r)**(-1/(alpha-1)) # Equation D.4 from Clauset, 2009
    g = np.sort(power)
    
    metric = 'KS'
   
    n_seeds = 10
    max_iterations = 1000


    alpha_history, x_min_history, cost_history = powerlaw_fit_SGD(g,max_iterations,learning_rate,n_seeds, momentum, fix_x_min=1, metric=metric, multiprocessing=True, early_stopping=False)
    alpha_history, x_min_history, cost_history = powerlaw_fit_SGD(g,max_iterations,learning_rate,n_seeds, momentum, fix_x_min=1, metric=metric, multiprocessing=False, early_stopping=False)
    # with np.printoptions(threshold=np.inf):
    #   print(np.asarray(alpha_history))


###############################################################################
# Stochastic DFN generation
###############################################################################

# def draw_lines2(image, lines):
#     for i in range(np.shape(lines)[0]):
#         image = cv2.line(image, (int(lines[i][0]), int(lines[i][1])),
#                         (int(lines[i][2]), int(lines[i][3])), (0, 0, 0), 1)
#     return image

def line_midpoint(line):
    x0, y0, x1, y1 = line
    y = (y0 + y1)/2
    x = (x0 + x1)/2
    return x, y

# Define a function to generate random samples from a power-law distribution
def powerlaw_sample(n, alpha, xmin):
    # Generate n uniform random numbers between 0 and 1
    u = np.random.uniform(0, 1, n)
    # Apply the inverse cumulative distribution function of the power-law distribution
    x = xmin * (1 - u) ** (-1 / (alpha - 1))
    return x

def fracture_center_generator(sizeX, sizeY, set_size):
    n = set_size
    positionX = np.random.random(n)
    positionY = np.random.random(n)

    positionX *= sizeX
    positionY *= sizeY

    return positionX, positionY, n


def fracture_center_to_line(center, theta, length):
    # print(np.shape(center), np.shape(length))
    X, Y = center

    # x0 = X + [-i * math.cos(math.radians(-90+theta))*0.5 for i in length]
    # y0 = Y + [-i * math.sin(math.radians(-90+theta))*0.5 for i in length]
    # x1 = X + [i * math.cos(math.radians(-90+theta))*0.5 for i in length]
    # y1 = Y + [i * math.sin(math.radians(-90+theta))*0.5 for i in length]

    x0 = X - length * np.cos(np.radians(-90+theta))*0.5
    y0 = Y - length * np.sin(np.radians(-90+theta))*0.5
    x1 = X + length * np.cos(np.radians(-90+theta))*0.5
    y1 = Y + length * np.sin(np.radians(-90+theta))*0.5
    return np.asarray([x0, y0, x1, y1]).T

def fracture_generator(set_size, theta, kappa, xmin, sizeX, sizeY, alpha, N):
    positionX, positionY, n = fracture_center_generator(sizeX,
                                                        sizeY, set_size)

    length = np.zeros((N, n))
    for i in range(0, N):
        # aux = atribute_generator(n, xmin, xmax, False, alpha = alpha)
        aux = powerlaw_sample(n, alpha, xmin)
        aux.sort()
        length[i] = aux

    theta_rad = np.radians(np.asarray([theta]*set_size))
    theta_deg = np.rad2deg(np.random.vonmises(theta_rad, kappa))
    mean_length = np.mean(length, axis=0)
    new_set = fracture_center_to_line((positionX, positionY), theta_deg, mean_length)
    return new_set

def generate_fractures(theta=0, alpha1=2.9, xmin=6.8, max_p21=1, dimensions=(500, 500), kappa = 100, monte_carlo_iterations=1000):
    sizeX, sizeY = dimensions
    count = 1
    while(True) :
        try:
            new_set1 = fracture_generator(count, theta, kappa, xmin, sizeX, sizeY, alpha1, monte_carlo_iterations)
            segm_group_angles = compute_line_angles(new_set1)
            sum_distances = np.sum(segm_group_angles[:,1])
            p21 = sum_distances/(sizeX*sizeY)
            if p21 <= max_p21:
                new_set0 = np.copy(new_set1)
                count += 1
                continue
            else:
                new_set1 = np.copy(new_set0)
                break
        except Exception as e:
            print(e)
        count += 1
    return new_set1
