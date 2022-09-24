import math
import seaborn as sns
import pandas
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

MAX_ITERATIONS = 1000
CHANGE_CUTOFF = 1 * (10 ** -20)
PH_CURVE_RESOLUTION = 0.01

Base = "Calcium Carbonate"
Salt = "Calcium Chloride"
Acid = "Hydrochloric acid"
OH = "Hydroxide"
H3O = "Hydronium"
Ions = "Calcium ions"

ACID = "Hydrochloric Acid dissociation"
BASE = "Tums dissociation"
ACID_BASE = "Acid Base Reaction"
WATER = "Water auto-ionization"

KW = 10.0 ** -8
K_ACID = 0.1 * (10.0 ** -2)
K_ACID_BASE = K_ACID / KW 
K_BASE = 0.45

REACTIONS = {
    ACID: (lambda chemicals: (- chemicals[H3O] - chemicals[Salt] - K_ACID + math.sqrt(
        ((chemicals[H3O] + chemicals[Salt] + K_ACID) ** 2) - 4 * (
                chemicals[H3O] * chemicals[Salt] - K_ACID * chemicals[Acid]))) / 2,
           [(Acid, False), (H3O, True), (Salt, True)]),
    BASE: (lambda chemicals: (- chemicals[OH] - chemicals[Ions] - K_BASE + math.sqrt(
        ((chemicals[OH] + chemicals[Ions] + K_BASE) ** 2) - 4 * (
                chemicals[Ions] * chemicals[OH] - K_BASE * chemicals[Base]))) / 2,
           [(Base, False), (OH, True), (Ions, True)]),
    WATER: (lambda chemicals: (- chemicals[OH] - chemicals[H3O] + math.sqrt(
        (chemicals[OH] + chemicals[H3O]) ** 2 - 4 * chemicals[OH] * chemicals[H3O] + 4 * KW)) / 2,
            [(OH, True), (H3O, True)]),
    ACID_BASE: (lambda chemicals: 2 / (- chemicals[OH] - chemicals[Acid] - K_ACID_BASE + math.sqrt(
        ((chemicals[OH] + chemicals[Acid] + K_ACID_BASE) ** 2) - 4 * (
                chemicals[OH] * chemicals[Acid] - K_ACID_BASE * chemicals[Salt]))),
                [(Acid, False), (OH, False), (Salt, True)])
}

INITIAL_CONCENTRATIONS = {
    OH: 10.0 ** (-7),
    Salt: 0,
    Acid: 1,
    H3O: 10.0 ** (-7),
    Ions: 0
}

sns.set()


class LinearGraph:
    def __init__(self, x_axis_name, y_axis_name, series_name):
        self.x_axis_name = x_axis_name
        self.y_axis_name = y_axis_name
        self.series_name = series_name

        self.x_values = []
        self.y_values = []
        self.series = []

        self.database = None

    def add_data_point(self, x_value, y_value, series):
        self.x_values.append(x_value)
        self.y_values.append(y_value)
        self.series.append(series)

    def add_data_for_each_series(self, x_value, series_with_data):
        for (series, value) in series_with_data.items():
            self.add_data_point(x_value, value, series)

    def _create_database(self):
        self.database = pandas.DataFrame(
            {self.x_axis_name: self.x_values, self.y_axis_name: self.y_values, self.series_name: self.series})

    def find_y_value(self, x_value, series):
        if self.database is None:
            self._create_database()

        return self.database.loc[(self.database[self.x_axis_name] == x_value) & (
                self.database[self.series_name] == series), self.y_axis_name].values[0]

    def graph(self, logarithmic_y_axis=False, title=None):
        if self.database is None:
            self._create_database()

        plot = sns.lineplot(x=self.x_axis_name, y=self.y_axis_name, hue=self.series_name, data=self.database)

        if logarithmic_y_axis:
            plt.yscale('symlog', linthreshy=10 ** -16)

        if title is not None:
            plot.set_title(title)

        plt.show()


def calculate_ph(hydronium_concentration):
    return -math.log10(hydronium_concentration)


def calculate_equilibrium(initial_naoh_concentration):
    concentrations = INITIAL_CONCENTRATIONS.copy()
    concentrations[Base] = initial_naoh_concentration

    total_changes = {}

    for name in REACTIONS.keys():
        total_changes[name] = 0

    for i in range(1, MAX_ITERATIONS):
        change_is_minimal = True

        for (name, (equation, changes)) in REACTIONS.items():
            required_change = equation(concentrations)
            if required_change != 0:
                for (chemical, positive) in changes:
                    if positive:
                        if required_change < 0:
                            required_change = - min(-required_change, concentrations[chemical])
                    else:
                        if required_change > 0:
                            required_change = min(required_change, concentrations[chemical])

                required_change /= 2

                for (chemical, positive) in changes:
                    if positive:
                        concentrations[chemical] += required_change
                    else:
                        concentrations[chemical] -= required_change

                total_changes[name] += required_change

                if required_change > CHANGE_CUTOFF:
                    change_is_minimal = False


    return concentrations, total_changes


def main():
    ph_graph = LinearGraph("Initial concentration of Base in Tums (mol/L)", "pH", "")
    change_graph = LinearGraph("Concentration of Base in Tums added (mol/L)", "Change in Concentration (mol/L)", "Reaction")

    for c in np.arange(0, 2, PH_CURVE_RESOLUTION):
        (chemicals, total_changes) = calculate_equilibrium(c)
        ph_graph.add_data_for_each_series(c, {"pH": calculate_ph(chemicals[H3O])})
        change_graph.add_data_for_each_series(c, total_changes)

    plt.ylim(0, 14)
    ph_graph.graph(title="pH as Tums is added to 0.1M Hydrochloric Acid")

    change_graph.graph(logarithmic_y_axis=False)


if __name__ == "__main__":
    main()