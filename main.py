import tkinter as tk
import itertools
from tkinter import messagebox
import numpy as np


class SimplexUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Simplex Solver")

        self.menubar = tk.Menu(root)
        help_menu = tk.Menu(self.menubar, tearoff=0)
        help_menu.add_command(label="Instrukcja", command=self.show_help)
        self.menubar.add_cascade(label="Pomoc", menu=help_menu)
        root.config(menu=self.menubar)

        header = tk.Label(root, text="Dane wyrobów", font=("Arial", 14, "bold"))
        header.grid(row=0, column=0, columnspan=5, pady=10)

        self.frame = tk.Frame(root)
        self.frame.grid(row=1, column=0, columnspan=5, padx=10)

        self.resource_names = []
        self.resource_headers = []
        self.resources = []
        self.rows = []

        self.product_header_label = tk.Label(self.frame, text="Wyrób", font=("Arial", 10, "bold"), width=20,
                                             anchor="center")
        self.product_header_label.grid(row=0, column=0, padx=3, pady=3)

        self.profit_header_label = tk.Label(self.frame, text="Zysk jednostkowy (c)", font=("Arial", 10, "bold"),
                                            width=20, anchor="center")

        self.initial_data = [
            (1, 1, 1, 17),
            (2, 3, 2, 43),
            (3, 1, 2, 20),
            (4, 1, 4, 24)
        ]

        btn_frame = tk.Frame(root)
        btn_frame.grid(row=2, column=0, columnspan=5, pady=10)

        tk.Button(btn_frame, text="➕ Dodaj wyrób", command=self.add_row).grid(row=0, column=0, padx=5)
        tk.Button(btn_frame, text="➖ Usuń ostatni", command=self.remove_row).grid(row=0, column=1, padx=5)

        resource_btn_frame = tk.Frame(root)
        resource_btn_frame.grid(row=3, column=0, columnspan=5, pady=5)

        tk.Button(resource_btn_frame, text="➕ Dodaj surowiec", command=self.add_resource).grid(row=0, column=0, padx=5)
        tk.Button(resource_btn_frame, text="➖ Usuń surowiec", command=self.remove_resource).grid(row=0, column=1,
                                                                                                 padx=5)

        calc_options_frame = tk.Frame(root)
        calc_options_frame.grid(row=4, column=0, columnspan=5, pady=10)

        opt_frame = tk.Frame(calc_options_frame)
        opt_frame.grid(row=0, column=0, padx=20)

        tk.Label(opt_frame, text="Optymalizacja:", font=("Arial", 10)).grid(row=0, column=0, padx=10)

        self.opt_type = tk.StringVar(value="max")
        tk.Radiobutton(opt_frame, text="Max", variable=self.opt_type, value="max").grid(row=0, column=1, padx=5, pady=2)
        tk.Radiobutton(opt_frame, text="Min", variable=self.opt_type, value="min").grid(row=0, column=2, padx=5, pady=2)

        combo_frame = tk.Frame(calc_options_frame)
        combo_frame.grid(row=0, column=1, padx=20)

        tk.Label(combo_frame, text="Wybierz ile wyrobów:").grid(row=0, column=0, padx=5)
        self.entry_combinations = tk.Entry(combo_frame, width=5)
        self.entry_combinations.grid(row=0, column=1, padx=5)
        self.entry_combinations.insert(0, "0")
        tk.Label(combo_frame, text="(0 = wszystkie)").grid(row=0, column=2, padx=5)

        tk.Label(root, text="Ograniczenia surowcowe (limity)", font=("Arial", 10, "bold")).grid(row=5, column=0,
                                                                                                columnspan=5,
                                                                                                pady=(10, 0))
        self.limits_frame = tk.Frame(root)
        self.limits_frame.grid(row=6, column=0, columnspan=5, pady=10)

        tk.Button(root, text="Oblicz", command=self.calculate).grid(row=7, column=0, columnspan=5, pady=15)

        self.add_resource(default_name="s1", default_val=30, default_type="<=")
        self.add_resource(default_name="s2", default_val=40, default_type="<=")

        for data in self.initial_data:
            self.add_row(default=data)

    def show_help(self):
        instructions = """
Instrukcja

1. Definiowanie Surowców:
    - Użyj przycisków ,,Dodaj surowiec" i ,,Usuń surowiec", aby zdefiniować liczbę surowców.
    - Wpisz maksymalną dostępną ilość każdego surowca w polach ,,Maks. zasób s1" i ,,Maks. zasób s2".

2. Definiowanie Wyrobów:
    - Użyj przycisków ,,Dodaj wyrób" i ,,Usuń wyrób", aby zdefiniować, ile wyrobów (zmiennych decyzyjnych) jest w problemie.
    - Dla każdego wyrobu, wpisz zużycie każdego surowca na jednostkę produktu.
    - W ostatniej kolumnie (,,Zysk jednostkowy") wpisz zysk (lub koszt) dla każdego wyrobu.

3. Typ Optymalizacji:
    - Wybierz ,,Max" dla maksymalizacji zysku lub ,,Min" dla minimalizacji kosztów.

4. Wybór Wyrobów:
    - Wpisz ,,0" (domyślnie), aby rozwiązać problem dla wszystkich wyrobów na liście.
    - Wpisz liczbę, aby przetestować każdą kombinację tej liczby wyrobów i znaleźć najlepszą z nich.

5. Obliczenia:
    - Naciśnij ,,Oblicz".

--- Ograniczenia Algorytmu ---
Obecna wersja algorytmu Simplex zakłada:
- Wszystkie ograniczenia są typu ,,mniejsze lub równe" (<=).
- Wszystkie limity surowców (prawa strona) są nieujemne (>= 0).
- Problemy z ograniczeniami ,,większe lub równe" (>=) lub ,,równe" (=) nie będą rozwiązane poprawnie i mogą prowadzić do błędnych wyników lub błędu aplikacji.
- Program wykrywa ,,rozwiązania nieograniczone", ale nie wykrywa ,,rozwiązań niedopuszczalnych", czyli sprzecznych ograniczeń.
        """
        messagebox.showinfo("Instrukcja", instructions)

    def add_resource(self, default_name=None, default_val=0, default_type="<="):
        res_idx = len(self.resource_names)
        resource_name = default_name if default_name else f"s{res_idx + 1}"
        col_idx = res_idx * 3

        lbl = tk.Label(self.limits_frame, text=f"{resource_name}:")
        lbl.grid(row=0, column=col_idx, padx=(10, 0), sticky='e')

        constraint_type_var = tk.StringVar(value=default_type)
        constraint_menu = tk.OptionMenu(self.limits_frame, constraint_type_var, "<=", ">=", "==")
        constraint_menu.config(width=3)
        constraint_menu.grid(row=0, column=col_idx + 1, padx=(5, 0))

        entry = tk.Entry(self.limits_frame, width=10)
        entry.insert(0, str(default_val))
        entry.grid(row=0, column=col_idx + 2, padx=(0, 10))

        self.resources.append({
            'label': lbl,
            'menu': constraint_menu,
            'constraint_type': constraint_type_var,
            'entry': entry
        })

        self.resource_names.append(resource_name)

        res_col_idx = res_idx + 1
        header_lbl = tk.Label(self.frame, text=resource_name, font=("Arial", 10, "bold"), width=20, anchor="center")
        header_lbl.grid(row=0, column=res_col_idx, padx=3, pady=3)
        self.resource_headers.append(header_lbl)

        self.profit_header_label.grid(row=0, column=res_col_idx + 1, padx=3, pady=3)

        for row_idx, row_widgets in enumerate(self.rows, start=1):
            e = tk.Entry(self.frame, width=20, justify="center")
            e.insert(0, "0")
            e.grid(row=row_idx, column=res_col_idx, padx=3, pady=3)
            row_widgets.insert(-1, e)
            profit_entry = row_widgets[-1]
            profit_entry.grid(row=row_idx, column=res_col_idx + 1, padx=3, pady=3)

    def remove_resource(self):
        if not self.resource_names:
            messagebox.showwarning("Błąd", "Brak surowców do usunięcia.")
            return

        r = self.resources.pop()
        r['label'].destroy()
        r['menu'].destroy()
        r['entry'].destroy()

        res_idx = len(self.resource_names)
        col_idx = res_idx + 1

        header_lbl = self.resource_headers.pop()
        header_lbl.destroy()

        self.profit_header_label.grid(row=0, column=col_idx, padx=3, pady=3)

        for row_idx, row_widgets in enumerate(self.rows, start=1):
            entry_to_remove = row_widgets.pop(-2)
            entry_to_remove.destroy()
            profit_entry = row_widgets[-1]
            profit_entry.grid(row=row_idx, column=col_idx, padx=3, pady=3)

        self.resource_names.pop()

    def add_row(self, default=None):
        idx = len(self.rows) + 1
        row_entries = []

        lbl = tk.Label(self.frame, text=f"Wyrób {idx}", width=20)
        lbl.grid(row=idx, column=0, padx=3, pady=3)
        row_entries.append(lbl)

        num_resources = len(self.resource_names)

        if default:
            values = default[1: 1 + num_resources]
            profit_val = default[1 + num_resources]
        else:
            values = [0] * num_resources
            profit_val = 0

        for col_idx, val in enumerate(values, start=1):
            e = tk.Entry(self.frame, width=20, justify="center")
            e.insert(0, str(val))
            e.grid(row=idx, column=col_idx, padx=3, pady=3)
            row_entries.append(e)

        e_profit = tk.Entry(self.frame, width=20, justify="center")
        e_profit.insert(0, str(profit_val))
        e_profit.grid(row=idx, column=num_resources + 1, padx=3, pady=3)
        row_entries.append(e_profit)

        self.rows.append(row_entries)

    def remove_row(self):
        if not self.rows:
            messagebox.showwarning("Błąd", "Brak wyrobów do usunięcia.")
            return
        for widget in self.rows.pop():
            widget.destroy()

    def _solve_simplex(self, products_data, b_vector_in, constraint_types_in, num_resources, opt_type):
        product_names = [p['name'] for p in products_data]
        c_vector = np.array([p['c'] for p in products_data], dtype=float)
        A_matrix = np.array([p['A_col'] for p in products_data], dtype=float).T

        b_vector = np.array(b_vector_in, dtype=float)
        constraint_types = list(constraint_types_in)

        num_decision_vars = len(products_data)

        for i in range(num_resources):
            if b_vector[i] < 0:
                b_vector[i] *= -1
                A_matrix[i, :] *= -1
                if constraint_types[i] == "<=":
                    constraint_types[i] = ">="
                elif constraint_types[i] == ">=":
                    constraint_types[i] = "<="

        slack_vars_count = 0
        surplus_vars_count = 0
        artificial_vars_count = 0

        slack_surplus_matrix = np.zeros((num_resources, num_resources))
        artificial_matrix = np.zeros((num_resources, num_resources))

        basis = []
        artificial_var_indices = []

        s_idx = 0
        a_idx = 0

        for i in range(num_resources):
            if constraint_types[i] == "<=":
                slack_surplus_matrix[i, s_idx] = 1
                basis.append(num_decision_vars + s_idx)
                s_idx += 1
                slack_vars_count += 1
            elif constraint_types[i] == ">=":
                slack_surplus_matrix[i, s_idx] = -1
                artificial_matrix[i, a_idx] = 1
                basis.append(num_decision_vars + num_resources + a_idx)
                artificial_var_indices.append(num_decision_vars + num_resources + a_idx)
                s_idx += 1
                a_idx += 1
                surplus_vars_count += 1
                artificial_vars_count += 1
            elif constraint_types[i] == "==":
                artificial_matrix[i, a_idx] = 1
                basis.append(num_decision_vars + num_resources + a_idx)
                artificial_var_indices.append(num_decision_vars + num_resources + a_idx)
                a_idx += 1
                artificial_vars_count += 1

        num_slack_surplus_vars = s_idx

        slack_surplus_matrix = slack_surplus_matrix[:, :num_slack_surplus_vars]
        artificial_matrix = artificial_matrix[:, :artificial_vars_count]

        tableau = np.hstack((A_matrix, slack_surplus_matrix, artificial_matrix))
        rhs = b_vector

        if artificial_vars_count > 0:
            c_phase1 = np.zeros(tableau.shape[1])
            c_phase1[num_decision_vars + num_slack_surplus_vars:] = 1.0

            cB = np.array([c_phase1[i] for i in basis])

            while True:
                z_j = cB @ tableau
                c_minus_z = c_phase1 - z_j

                if np.all(c_minus_z >= -1e-9):
                    break

                pivot_col_idx = np.argmin(c_minus_z)
                pivot_column = tableau[:, pivot_col_idx]

                if np.all(pivot_column <= 1e-9):
                    return {"status": "unbounded", "profit": -np.inf, "product_names": product_names}

                ratios = np.full(num_resources, np.inf)
                for i in range(num_resources):
                    if pivot_column[i] > 1e-9:
                        ratios[i] = rhs[i] / pivot_column[i]

                pivot_row_idx = np.argmin(ratios)
                pivot_element = tableau[pivot_row_idx, pivot_col_idx]

                basis[pivot_row_idx] = pivot_col_idx
                cB[pivot_row_idx] = c_phase1[pivot_col_idx]

                pivot_row = tableau[pivot_row_idx, :]
                tableau[pivot_row_idx, :] = pivot_row / pivot_element
                rhs[pivot_row_idx] = rhs[pivot_row_idx] / pivot_element

                for i in range(num_resources):
                    if i != pivot_row_idx:
                        factor = tableau[i, pivot_col_idx]
                        tableau[i, :] -= factor * tableau[pivot_row_idx, :]
                        rhs[i] -= factor * rhs[pivot_row_idx]

            final_profit_phase1 = cB @ rhs
            if final_profit_phase1 > 1e-6:
                return {"status": "infeasible"}

            cols_to_delete = list(range(num_decision_vars + num_slack_surplus_vars, tableau.shape[1]))
            tableau = np.delete(tableau, cols_to_delete, axis=1)

        c_phase2 = np.concatenate((c_vector, np.zeros(num_slack_surplus_vars)))

        if opt_type == "min":
            c_phase2 = -c_phase2

        cB = np.array([c_phase2[i] for i in basis])
        c_minus_z = np.zeros(tableau.shape[1])  # Inicjalizacja

        while True:
            z_j = cB @ tableau
            c_minus_z = c_phase2 - z_j

            if np.all(c_minus_z <= 1e-9):
                break

            pivot_col_idx = np.argmax(c_minus_z)
            pivot_column = tableau[:, pivot_col_idx]

            if np.all(pivot_column <= 1e-9):
                return {
                    "status": "unbounded",
                    "profit": np.inf if opt_type == "max" else -np.inf,
                    "product_names": product_names
                }

            ratios = np.full(num_resources, np.inf)
            for i in range(num_resources):
                if pivot_column[i] > 1e-9:
                    ratios[i] = rhs[i] / pivot_column[i]

            pivot_row_idx = np.argmin(ratios)
            if np.all(ratios == np.inf):
                return {
                    "status": "unbounded",
                    "profit": np.inf if opt_type == "max" else -np.inf,
                    "product_names": product_names
                }

            pivot_element = tableau[pivot_row_idx, pivot_col_idx]

            basis[pivot_row_idx] = pivot_col_idx
            cB[pivot_row_idx] = c_phase2[pivot_col_idx]

            pivot_row = tableau[pivot_row_idx, :]
            tableau[pivot_row_idx, :] = pivot_row / pivot_element
            rhs[pivot_row_idx] = rhs[pivot_row_idx] / pivot_element

            for i in range(num_resources):
                if i != pivot_row_idx:
                    factor = tableau[i, pivot_col_idx]
                    tableau[i, :] -= factor * tableau[pivot_row_idx, :]
                    rhs[i] -= factor * rhs[pivot_row_idx]

        all_vars_indices = set(range(tableau.shape[1]))
        non_basic_indices = all_vars_indices - set(basis)

        has_alternatives = False
        for idx in non_basic_indices:
            if abs(c_minus_z[idx]) < 1e-9:
                if np.any(tableau[:, idx] > 1e-9):
                    has_alternatives = True
                    break

        final_profit = cB @ rhs
        if opt_type == "min":
            final_profit = -final_profit

        solution_values = np.zeros(num_decision_vars, dtype=float)
        for i, var_idx in enumerate(basis):
            if var_idx < num_decision_vars:
                solution_values[var_idx] = rhs[i]

        return {
            "status": "optimal",
            "profit": final_profit,
            "product_names": product_names,
            "solution_values": solution_values,
            "has_alternatives": has_alternatives
        }

    def _display_result(self, result, opt_type, k_val):

        if result['status'] == 'unbounded':
            product_names_str = ", ".join(result['product_names'])
            messagebox.showwarning("Rozwiązanie nieograniczone",
                                   f"Problem nie ma ograniczonego rozwiązania optymalnego.\n"
                                   f"(Dotyczy kombinacji: {product_names_str})")
            return

        if result['status'] == 'infeasible':
            messagebox.showerror("Brak rozwiązania",
                                 "Rozwiązanie niedopuszczalne!\n\n"
                                 "Ograniczenia są ze sobą sprzeczne (np. x <= 5 i x >= 10).")
            return

        if result['status'] == 'optimal':
            title = f"Wynik optymalny ({opt_type})"
            profit_str = f"{result['profit']:.2f} zł"

            solution_str = ""
            for name, val in zip(result['product_names'], result['solution_values']):
                solution_str += f"{name}: {val:.2f} jedn.\n"

            if k_val > 0:
                product_names_str = ", ".join(result['product_names'])
                header = f"Najlepsza kombinacja {k_val} wyrobów: ({product_names_str})\n\n"
            else:
                header = "Osiągnięto optymalne rozwiązanie!\n\n"

            message = header + f"Wartość funkcji celu (Zysk/Koszt): {profit_str}\n\n" \
                               f"Należy produkować:\n{solution_str}"

            if result.get("has_alternatives", False):
                message += "\n\nUwaga: Istnieją alternatywne rozwiązania optymalne (rozwiązaniem jest odcinek), które dają ten sam wynik."

            messagebox.showinfo(title, message)

        elif result['status'] == 'error':
            messagebox.showerror("Błąd", result['message'])

    def calculate(self):
        try:
            b_vector = [float(r['entry'].get()) for r in self.resources]
            constraint_types = [r['constraint_type'].get() for r in self.resources]
            num_resources = len(self.resources)
            opt_type = self.opt_type.get()

            if len(self.rows) < 1:
                messagebox.showerror("Błąd", "Potrzebny jest co najmniej 1 wyrób.")
                return

            all_products_data = []
            for row_widgets in self.rows:
                name = row_widgets[0]["text"]
                c = float(row_widgets[-1].get())

                resource_cols = []
                for i in range(num_resources):
                    resource_cols.append(float(row_widgets[i + 1].get()))

                all_products_data.append({"name": name, "c": c, "A_col": resource_cols})

            k_str = self.entry_combinations.get()
            k = int(k_str) if k_str.isdigit() else 0

            results_list = []

            if k > 0:
                if k > len(all_products_data):
                    messagebox.showerror("Błąd", f"Nie można wybrać {k} wyrobów z {len(all_products_data)} dostępnych.")
                    return

                if not all_products_data:
                    messagebox.showerror("Błąd", "Brak wyrobów do przetworzenia.")
                    return

                for product_combo in itertools.combinations(all_products_data, k):
                    result = self._solve_simplex(list(product_combo), b_vector, constraint_types, num_resources,
                                                 opt_type)
                    results_list.append(result)

                if not results_list:
                    messagebox.showinfo("Brak wyników", "Nie udało się obliczyć żadnej kombinacji.")
                    return

                valid_results = [r for r in results_list if r['status'] in ('optimal', 'unbounded')]

                if not valid_results:
                    if all(r['status'] == 'infeasible' for r in results_list):
                        messagebox.showerror("Błąd",
                                             "Wszystkie możliwe kombinacje produktów prowadzą do sprzecznych ograniczeń (brak rozwiązania).")
                    else:
                        messagebox.showerror("Błąd", "Wszystkie kombinacje zakończyły się błędem.")
                    return

                if opt_type == "max":
                    best_result = max(valid_results, key=lambda x: x.get("profit", -np.inf))
                else:
                    best_result = min(valid_results, key=lambda x: x.get("profit", np.inf))

                self._display_result(best_result, opt_type, k)

            else:
                if not all_products_data:
                    messagebox.showerror("Błąd", "Brak wyrobów do przetworzenia.")
                    return

                result = self._solve_simplex(all_products_data, b_vector, constraint_types, num_resources, opt_type)
                self._display_result(result, opt_type, 0)

        except ValueError:
            messagebox.showerror("Błąd danych",
                                 "Wprowadzono niepoprawne dane. Upewnij się, że wszystkie pola zawierają liczby.")
        except Exception as e:
            messagebox.showerror("Błąd obliczeń", f"Wystąpił nieoczekiwany błąd: {e}")


if __name__ == "__main__":
    root = tk.Tk()
    app = SimplexUI(root)
    root.mainloop()