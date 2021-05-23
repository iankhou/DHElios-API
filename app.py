import time
from flask import Flask, Response
import io
import random
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import *
import numpy as np
from numpy import *
from sympy import *
from sympy import symbols
import math


app = Flask(__name__)


@app.route('/')
def welcome():
    return 'Welcome to DHElios'


@app.route('/time')
def get_current_time():
    return {'time': time.time()}


@app.route('/graph')
def get_graph():
    return None


@app.route('/plot.png')
def plot_png():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')


lowest_max_power = 10
highest_max_power = 120


def create_max_power_list():
    max_power_list = []
    i = lowest_max_power
    while lowest_max_power <= i and i <= highest_max_power:
        max_power_list.append(i)
        i = i + 40
    return max_power_list


max_power_list = create_max_power_list()

# constants
# inputed panel areas list
panel_num_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# this is the number of seconds that we are evaluating at for delta T (12 hours)
evaluate_num = 3600*12  # seconds in hour times number of hours

# resulting panel area list - each of these will be the "x's" in the function
single_panel_area = 2
panel_area_list = []
for item in range(0, len(panel_num_list)):
    new_panel_area = panel_num_list[item] * single_panel_area
    panel_area_list.append(new_panel_area)

# # initial variable values
slope = -25.1  # slope
T_air = 303  # air temperature, units: K
y_intercept = 4100  # y-intercept
k = 0.022  # conductivity of polyurethane, units: W/(m*K)
surface_area = 6  # of water tank, units: m^2
thickness = 0.03  # depth/thickness of polyurethane, untis: m
initial_water_temp = 300  # initial water temp, units: K

# define constants for mass and heat capacity
mass = 1000  # kg for 10k
heat_capacity = 4179.6  # check units #

# THE X OF THE FUNCTION IS CURRENT WATER TEMP, WHICH THEN BECOMES TIME
x = symbols('x')

coil_efficiency = 1
max_power = 1000

solar_output = 1760


# TOTAL COST OF PV PANEL SYSTEM
# define constants for cost
PV_fixed_costs = 4000  # need a real estimate here - these are fake numbers
PV_panel_cost = 600  # also need a real estimate

# ANNUAL SAVINGS
total_firewood_cost = 6000  # in dollars
fraction_wood_used = 2/3    # fraction of wood that can be adjusted
final_water_temp = 373      # in Kelvin

# TOTAL COST OF GLAZED PLATE SYSTEM
# define constants for cost
fixed_costs = 3000  # need a real estimate here - these are fake numbers
panel_cost = 400  # also need a real estimate


# ANNUAL SAVINGS
total_firewood_cost = 6000  # in dollars
fraction_wood_used = 2/3    # fraction of wood that can be adjusted
final_water_temp = 373      # in Kelvin


# the x is time, y is temperature
def glazed_plate_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area) -> float:
    a = (slope*panel_area) - ((k*surface_area)/thickness)
    b = (slope*-1*T_air*panel_area) + (y_intercept*panel_area) + \
        (((k*surface_area)/thickness)*T_air)
    constant = 298*a + b
    time_expression = (1/a)*(constant*e**((x*a)/(mass*heat_capacity)) - b)
    return time_expression


def cooling_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area) -> float:
    glazed_plate_evaluation = glazed_plate_function(
        evaluate_num, T_air, k, surface_area, thickness, initial_water_temp, panel_area)
    constant = glazed_plate_evaluation - T_air
    a = ((k*surface_area)/thickness)/(mass*heat_capacity)
    time_expression = (constant * (e**(-1*((x-evaluate_num)*a)))) + T_air
    return time_expression


def glazed_piecewise(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area):
    glazed_time_expression = Piecewise((glazed_plate_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area),
                                       x <= evaluate_num), (cooling_function(x, T_air, k, surface_area, thickness, initial_water_temp, panel_area), x > evaluate_num))
    return glazed_time_expression


def PV_panel_function(x, surface_area, thickness, solar_output, max_power) -> float:
    a = (k*surface_area)/thickness
    b = max_power * solar_output * coil_efficiency
    constant = (298-T_air)*a - b
    PV_time_expression = T_air + (constant * e**(-a*x) - b)/(a)
    return PV_time_expression


def PV_cooling_function(x, surface_area, thickness, solar_output, max_power) -> float:
    PV_evaluation = PV_panel_function(
        x, surface_area, thickness, solar_output, max_power)
    constant = PV_evaluation - T_air
    a = ((k*surface_area)/thickness)/(mass*heat_capacity)
    time_expression = (constant * (e**(-1*((x-evaluate_num)*a)))) + T_air
    return time_expression


def PV_piecewise(x, surface_area, thickness, solar_output, max_power):
    PV_time_expression = Piecewise((PV_panel_function(x, surface_area, thickness, solar_output, max_power), x <=
                                   evaluate_num), (PV_cooling_function(x, surface_area, thickness, solar_output, max_power), x > evaluate_num))
    return PV_time_expression

# loop through each of the panel areas and get a resulting list of expressions of time for each plate area


def new_time_expression_list(x, T_air, k, surface_area, thickness, initial_water_temp):
    time_expression_list = []
    for item in range(0, len(panel_num_list)):
        glazed_new_time_exp = glazed_piecewise(
            x, T_air, k, surface_area, thickness, initial_water_temp, panel_area_list[item])
        time_expression_list.append(glazed_new_time_exp)
    return time_expression_list


def PV_new_time_expression_list(x, surface_area, thickness, solar_output):
    time_expression_list = []
    for item in range(0, len(max_power_list)):
        PV_new_time_exp = PV_piecewise(
            x, surface_area, thickness, solar_output, max_power_list[item])
        time_expression_list.append(PV_new_time_exp)
    return time_expression_list

# calculate delta_T for each plate area


def new_delta_T_list(time_expression_list):
    delta_T_list = []
    at_final_time = new_time_expression_list(
        evaluate_num*2, T_air, k, surface_area, thickness, initial_water_temp)
    at_initial_time = new_time_expression_list(
        0, T_air, k, surface_area, thickness, initial_water_temp)
    for index in range(0, len(panel_num_list)):
        delta_T_list.append(at_final_time[index] - at_initial_time[index])
    return delta_T_list


def PV_new_delta_T_list(time_expression_list):
    PV_delta_T_list = []
    at_final_time = PV_new_time_expression_list(
        evaluate_num*2, surface_area, thickness, solar_output)
    at_initial_time = PV_new_time_expression_list(
        0, surface_area, thickness, solar_output)
    for index in range(0, len(max_power_list)):
        PV_delta_T_list.append(at_final_time[index] - at_initial_time[index])
    return PV_delta_T_list

# function to calculate cost of the system


def PV_calculate_system_cost(PV_panel_cost, max_power, PV_fixed_costs):
    PV_total_cost = PV_panel_cost * max_power + PV_fixed_costs
    return PV_total_cost

# loop through all of the panel numbers to get a list of total costs


def PV_new_total_cost_list():
    PV_total_cost_list = []
    for item in max_power_list:
        new_cost = PV_calculate_system_cost(
            PV_panel_cost, item, PV_fixed_costs)
        PV_total_cost_list.append(new_cost)
    return PV_total_cost_list

# function to calculate the annual savings


def PV_calculate_annual_savings(PV_delta_T, fraction_wood_used):
    annual_savings = fraction_wood_used * \
        total_firewood_cost * (PV_delta_T)/final_water_temp
    return annual_savings

# loop through all the delta Ts to get a list of annual savings for each delta T


def PV_new_annual_savings_list(PV_delta_T_list):
    PV_annual_savings_list = []
    for item in PV_delta_T_list:
        new_annual_savings = PV_calculate_annual_savings(
            item, fraction_wood_used)
        PV_annual_savings_list.append(new_annual_savings)

    return PV_annual_savings_list

# PAYBACK PERIOD
# function to calculate payback period


def PV_calculate_payback_period(total_cost, annual_savings):
    payback_period = total_cost/annual_savings
    return payback_period

# loop through all the total costs to get a list of payback periods for each panel number


def PV_new_payback_period_list(PV_total_cost_list, PV_annual_savings_list):
    PV_payback_period_list = []
    for item in range(0, len(max_power_list)):
        new_payback = PV_calculate_payback_period(
            PV_total_cost_list[item], PV_annual_savings_list[item])
        PV_payback_period_list.append(new_payback)
    return PV_payback_period_list

# function to calculate cost of the system


def calculate_system_cost(panel_cost, panel_num, fixed_costs):
    total_cost = panel_cost * panel_num + fixed_costs
    return total_cost

# loop through all of the panel numbers to get a list of total costs


def new_total_cost_list():
    total_cost_list = []
    for item in panel_num_list:
        new_cost = calculate_system_cost(panel_cost, item, fixed_costs)
        total_cost_list.append(new_cost)

    return total_cost_list

# function to calculate the annual savings


def calculate_annual_savings(delta_T, fraction_wood_used):
    annual_savings = fraction_wood_used * \
        total_firewood_cost * (delta_T)/final_water_temp
    return annual_savings

# loop through all the delta Ts to get a list of annual savings for each delta T


def new_annual_savings_list(delta_T_list):
    annual_savings_list = []
    for item in delta_T_list:
        new_annual_savings = calculate_annual_savings(item, fraction_wood_used)
        annual_savings_list.append(new_annual_savings)

    return annual_savings_list

# PAYBACK PERIOD
# function to calculate payback period


def calculate_payback_period(total_cost, annual_savings):
    payback_period = total_cost/annual_savings
    return payback_period

# loop through all the total costs to get a list of payback periods for each panel number


def new_payback_period_list(total_cost_list, annual_savings_list):
    payback_period_list = []
    for item in range(0, len(panel_num_list)):
        new_payback = calculate_payback_period(
            total_cost_list[item], annual_savings_list[item])
        payback_period_list.append(new_payback)
    return payback_period_list


def plot_time_temp_graph(T_air, k, surface_area, thickness, initial_water_temp, fig):
    # axes for temp vs time graph
    sub1_x_axis_0 = 0
    sub1_x_axis_f = evaluate_num * 2  # this should be 12 hours, in hours
    sub1_y_axis_0 = 280
    sub1_y_axis_f = 500  # boiling temp in K

    sub1 = fig.add_subplot(2, 2, 1)
    sub1.cla()  # clear the subplot's axes before drawing on it
    sub1.set_title('glazed temp vs. time')
    sub1.set_xlabel('time (sec)')
    sub1.set_ylabel('temperature (K)')
    x = np.linspace(0, evaluate_num * 2, 200)

    sub1.set_xlim([sub1_x_axis_0, sub1_x_axis_f])
    sub1.set_ylim([sub1_y_axis_0, sub1_y_axis_f])

    temp_time_graph_list = []

    for item in range(0, len(panel_num_list)):
        pw = np.piecewise(x, [x < evaluate_num, x >= evaluate_num], [
            lambda z: glazed_plate_function(z, T_air, k, surface_area, thickness, initial_water_temp,
                                            panel_area_list[item]), lambda q: cooling_function(q, T_air, k, surface_area, thickness, initial_water_temp, panel_area_list[item])])
        temp_time_graph_list.append(sub1.plot(x, pw))


def PV_plot_time_temp_graph(surface_area, thickness, solar_output):
    # axes for temp vs time graph
    sub1_x_axis_0 = 0
    sub1_x_axis_f = evaluate_num * 2  # this should be 12 hours, in hours

    sub1 = plt.subplot(2, 2, 1)
    sub1.cla()  # clear the subplot's axes before drawing on it
    sub1.set_title('PV temp vs. time')
    sub1.set_xlabel('time (sec)')
    sub1.set_ylabel('temperature (K)')

    x = np.linspace(0, evaluate_num * 2, 200)

    sub1.set_xlim([sub1_x_axis_0, sub1_x_axis_f])

    temp_time_graph_list = []

    for item in range(0, len(max_power_list)):
        pw = np.piecewise(x, [x < evaluate_num, x >= evaluate_num], [
            lambda z: PV_panel_function(
                z, surface_area, thickness, solar_output, max_power_list[item]),
            lambda q: PV_cooling_function(q, surface_area, thickness, solar_output, max_power_list[item])])
        temp_time_graph_list.append(sub1.plot(x, pw))


def plot_payback_graph(delta_T_output, payback_period_output, fig):
    # axes for payback period vs delta T graph
    sub2_x_axis_0 = 2
    sub2_x_axis_f = 7     # hopefully not more than 10 years
    sub2_y_axis_0 = 0
    sub2_y_axis_f = 200  # boiling temp in K

    sub2 = fig.add_subplot(2, 2, 2)
    sub2.set_title('payback vs. delta T')
    sub2.set_yticks(np.arange(sub2_x_axis_0, sub2_x_axis_f, step=1))
    sub2.set_xticks(np.arange(sub2_y_axis_0, sub2_y_axis_f, step=20))
    sub2.set_ylabel('payback period (yrs)')
    sub2.set_xlabel('delta T (K)')
    sub2.set_ylim([sub2_x_axis_0, sub2_x_axis_f])
    sub2.set_xlim([sub2_y_axis_0, sub2_y_axis_f])
    f_d2, = sub2.plot(delta_T_output, payback_period_output, marker='o')


def PV_plot_payback_graph(PV_delta_T_output, PV_payback_period_output):
    # axes for payback period vs delta T graph
    sub2_x_axis_0 = 0
    sub2_x_axis_f = 10     # hopefully not more than 10 years
    sub2_y_axis_0 = 0
    sub2_y_axis_f = 373.16  # boiling temp in K

    sub2 = plt.subplot(2, 2, 2)
    sub2.set_title('PV payback vs. delta T')
    sub2.set_xticks(np.arange(sub2_x_axis_0, sub2_x_axis_f, step=5))
    sub2.set_yticks(np.arange(sub2_y_axis_0, sub2_y_axis_f, step=50))
    sub2.set_xlabel('payback period (yrs)')
    sub2.set_ylabel('delta T (K)')


def create_figure():
    time_expression_output = new_time_expression_list(
        x, T_air, k, surface_area, thickness, initial_water_temp)
    PV_time_expression_output = PV_new_time_expression_list(
        x, surface_area, thickness, solar_output)
    delta_T_output = new_delta_T_list(time_expression_output)
    PV_delta_T_output = PV_new_delta_T_list(PV_time_expression_output)
    PV_total_cost_output = PV_new_total_cost_list()
    PV_annual_savings_output = PV_new_annual_savings_list(PV_delta_T_output)
    PV_payback_period_output = PV_new_payback_period_list(
        PV_total_cost_output, PV_annual_savings_output)
    total_cost_output = new_total_cost_list()
    annual_savings_output = new_annual_savings_list(delta_T_output)
    payback_period_output = new_payback_period_list(
        total_cost_output, annual_savings_output)
    complete_array = np.array([panel_num_list, panel_area_list, delta_T_output,
                              total_cost_output, annual_savings_output, payback_period_output])
    # fig = plt.figure(figsize=(10, 6))
    fig = Figure()

    plot_time_temp_graph(T_air, k, surface_area,
                         thickness, initial_water_temp, fig)
    plot_payback_graph(delta_T_output, payback_period_output, fig)

    # fig.show()
    # FigureCanvas(fig).draw()
    return fig


create_figure()
