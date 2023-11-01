from .solar import (
    estimate_panel_cost, estimate_battery_cost, estimate_balance_cost)
from .mini_grid import estimate_lv_connection_cost
from . import (
    prepare_component_cost_by_year, prepare_external_cost,
    prepare_internal_cost)


def estimate_internal_cost(**keywords):
    return prepare_internal_cost([
        estimate_system_capacity_cost,
        estimate_internal_distribution_cost,
    ], keywords)


def estimate_external_cost(**keywords):
    return prepare_external_cost([
    ], keywords)

def estimate_internal_distribution_cost(**keywords):
    component_cost_by_year, d = prepare_component_cost_by_year([
        ('lv_connection', estimate_solar_home_lv_connection_cost),
    ], keywords, prefix='solar_home_')
    d['internal_distribution_cost_by_year'] = component_cost_by_year
    return d


def estimate_system_capacity_cost(**keywords):
    component_cost_by_year, d = prepare_component_cost_by_year([
        ('panel', estimate_solar_home_panel_cost),
        ('battery', estimate_solar_home_battery_cost),
        ('balance', estimate_solar_home_balance_cost),
    ], keywords, prefix='solar_home_')
    d['system_capacity_cost_by_year'] = component_cost_by_year
    return d


def estimate_solar_home_panel_cost(
        final_consumption_in_kwh_per_year,
        peak_hours_of_sun_per_year,
        solar_home_system_loss_as_percent_of_total_production,
        solar_home_panel_table):
    return estimate_panel_cost(
        final_consumption_in_kwh_per_year,
        peak_hours_of_sun_per_year,
        solar_home_system_loss_as_percent_of_total_production,
        solar_home_panel_table)


def estimate_solar_home_battery_cost(
        solar_home_panel_actual_system_capacity_in_kw,
        solar_home_battery_kwh_per_panel_kw,
        solar_home_battery_raw_cost_per_battery_kwh,
        solar_home_battery_installation_cost_as_percent_of_raw_cost,
        solar_home_battery_maintenance_cost_per_year_as_percent_of_raw_cost,
        solar_home_battery_lifetime_in_years):
    return estimate_battery_cost(
        solar_home_panel_actual_system_capacity_in_kw,
        solar_home_battery_kwh_per_panel_kw,
        solar_home_battery_raw_cost_per_battery_kwh,
        solar_home_battery_installation_cost_as_percent_of_raw_cost,
        solar_home_battery_maintenance_cost_per_year_as_percent_of_raw_cost,
        solar_home_battery_lifetime_in_years)


def estimate_solar_home_balance_cost(
        solar_home_panel_actual_system_capacity_in_kw,
        solar_home_balance_raw_cost_per_panel_kw,
        solar_home_balance_installation_cost_as_percent_of_raw_cost,
        solar_home_balance_maintenance_cost_per_year_as_percent_of_raw_cost,
        solar_home_balance_lifetime_in_years):
    return estimate_balance_cost(
        solar_home_panel_actual_system_capacity_in_kw,
        solar_home_balance_raw_cost_per_panel_kw,
        solar_home_balance_installation_cost_as_percent_of_raw_cost,
        solar_home_balance_maintenance_cost_per_year_as_percent_of_raw_cost,
        solar_home_balance_lifetime_in_years)

def estimate_solar_home_lv_connection_cost(
        final_connection_count,
        solar_home_lv_connection_raw_cost,
        solar_home_lv_connection_installation_cost_as_percent_of_raw_cost,
        solar_home_lv_connection_maintenance_cost_per_year_as_percent_of_raw_cost,  # noqa
        solar_home_lv_connection_lifetime_in_years):
    return estimate_lv_connection_cost(
        final_connection_count,
        solar_home_lv_connection_raw_cost,
        solar_home_lv_connection_installation_cost_as_percent_of_raw_cost,
        solar_home_lv_connection_maintenance_cost_per_year_as_percent_of_raw_cost,  # noqa
        solar_home_lv_connection_lifetime_in_years)

