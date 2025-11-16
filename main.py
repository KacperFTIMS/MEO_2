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

        tk.Button(btn_frame, text="Dodaj wyrób", command=self.add_row).grid(row=0, column=0, padx=5)
        tk.Button(btn_frame, text="Usuń wyrób", command=self.remove_row).grid(row=0, column=1, padx=5)

        resource_btn_frame = tk.Frame(root)
        resource_btn_frame.grid(row=3, column=0, columnspan=5, pady=5)

        tk.Button(resource_btn_frame, text="Dodaj surowiec", command=self.add_resource).grid(row=0, column=0, padx=5)
        tk.Button(resource_btn_frame, text="Usuń surowiec", command=self.remove_resource).grid(row=0, column=1,
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

        tk.Label(combo_frame, text="Wybierz wyrobów:").grid(row=0, column=0, padx=5)
        self.entry_combinations = tk.Entry(combo_frame, width=5)
        self.entry_combinations.grid(row=0, column=1, padx=5)
        self.entry_combinations.insert(0, "0")
        tk.Label(combo_frame, text="(0 = wszystkie)").grid(row=0, column=2, padx=5)

        self.limits_frame = tk.Frame(root)
        self.limits_frame.grid(row=5, column=0, columnspan=5, pady=10)

        tk.Button(root, text="Oblicz", command=self.calculate).grid(row=6, column=0, columnspan=5, pady=15)

        self.add_resource(default_name="s1", default_val=30)
        self.add_resource(default_name="s2", default_val=40)

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

    def add_resource(self, default_name=None, default_val=0):
        res_idx = len(self.resource_names)
        resource_name = default_name if default_name else f"s{res_idx + 1}"

        col_idx = res_idx + 1
        header_lbl = tk.Label(self.frame, text=resource_name, font=("Arial", 10, "bold"), width=20, anchor="center")
        header_lbl.grid(row=0, column=col_idx, padx=3, pady=3)
        self.resource_headers.append(header_lbl)

        self.profit_header_label.grid(row=0, column=col_idx + 1, padx=3, pady=3)

        lbl = tk.Label(self.limits_frame, text=f"Maks. zasób {resource_name}:")
        lbl.grid(row=0, column=2 * res_idx, padx=5)
        entry = tk.Entry(self.limits_frame, width=10)
        entry.insert(0, str(default_val))
        entry.grid(row=0, column=2 * res_idx + 1, padx=5)
        self.resources.append({'label': lbl, 'entry': entry})

        for row_idx, row_widgets in enumerate(self.rows, start=1):
            e = tk.Entry(self.frame, width=20, justify="center")
            e.insert(0, "0")
            e.grid(row=row_idx, column=col_idx, padx=3, pady=3)

            row_widgets.insert(-1, e)

            profit_entry = row_widgets[-1]
            profit_entry.grid(row=row_idx, column=col_idx + 1, padx=3, pady=3)

        self.resource_names.append(resource_name)

    def remove_resource(self):
        if not self.resource_names:
            messagebox.showwarning("Błąd", "Brak surowców do usunięcia.")
            return

        res_idx = len(self.resource_names) - 1
        col_idx = res_idx + 1

        header_lbl = self.resource_headers.pop()
        header_lbl.destroy()

        self.profit_header_label.grid(row=0, column=col_idx, padx=3, pady=3)

        r = self.resources.pop()
        r['label'].destroy()
        r['entry'].destroy()

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

    def _solve_simplex(self, products_data, b_vector, num_resources, opt_type):
        product_names = [p['name'] for p in products_data]
        c_vector = [p['c'] for p in products_data]
        A_matrix_T = [p['A_col'] for p in products_data]

        num_decision_vars = len(products_data)
        num_slack_vars = num_resources
        num_total_vars = num_decision_vars + num_slack_vars

        c_j_row = np.array(c_vector + [0] * num_slack_vars, dtype=float)

        if not A_matrix_T:
            return {"status": "error", "message": "Brak wyrobów do przetworzenia."}

        A_matrix = np.array(A_matrix_T).T
        I_matrix = np.eye(num_resources)
        constraints_matrix = np.hstack((A_matrix, I_matrix))

        rhs_b = np.array(b_vector, dtype=float)

        cB = np.zeros(num_resources, dtype=float)

        basic_var_indices = list(range(num_decision_vars, num_total_vars))

        while True:
            z_j = cB @ constraints_matrix
            c_j_minus_z_j = c_j_row - z_j

            if (opt_type == "max" and np.all(c_j_minus_z_j <= 0)) or \
                    (opt_type == "min" and np.all(c_j_minus_z_j >= 0)):
                break

            if opt_type == "max":
                pivot_col_idx = np.argmax(c_j_minus_z_j)
            else:
                pivot_col_idx = np.argmin(c_j_minus_z_j)

            pivot_column = constraints_matrix[:, pivot_col_idx]

            if np.all(pivot_column <= 0):
                return {
                    "status": "unbounded",
                    "profit": np.inf if opt_type == "max" else -np.inf,
                    "product_names": product_names
                }

            ratios = np.full(num_resources, np.inf)
            for i in range(num_resources):
                if pivot_column[i] > 0:
                    ratios[i] = rhs_b[i] / pivot_column[i]

            pivot_row_idx = np.argmin(ratios)
            pivot_element = constraints_matrix[pivot_row_idx, pivot_col_idx]

            cB[pivot_row_idx] = c_j_row[pivot_col_idx]
            basic_var_indices[pivot_row_idx] = pivot_col_idx

            pivot_row = constraints_matrix[pivot_row_idx, :]
            pivot_row /= pivot_element
            rhs_b[pivot_row_idx] /= pivot_element

            for i in range(num_resources):
                if i != pivot_row_idx:
                    factor = constraints_matrix[i, pivot_col_idx]
                    constraints_matrix[i, :] -= factor * pivot_row
                    rhs_b[i] -= factor * rhs_b[pivot_row_idx]

        final_profit = cB @ rhs_b

        solution_values = np.zeros(num_decision_vars, dtype=float)
        for i, var_idx in enumerate(basic_var_indices):
            if var_idx < num_decision_vars:  # Jeśli zmienna bazowa jest zmienną decyzyjną
                solution_values[var_idx] = rhs_b[i]

        return {
            "status": "optimal",
            "profit": final_profit,
            "product_names": product_names,
            "solution_values": solution_values
        }

    def _display_result(self, result, opt_type, k_val):

        if result['status'] == 'unbounded':
            product_names_str = ", ".join(result['product_names'])
            messagebox.showwarning("Rozwiązanie nieograniczone",
                                   f"Problem nie ma ograniczonego rozwiązania optymalnego.\n"
                                   f"(Dotyczy kombinacji: {product_names_str})")
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

            messagebox.showinfo(title, message)

        elif result['status'] == 'error':
            messagebox.showerror("Błąd", result['message'])

    def calculate(self):
        try:
            b_vector = [float(r['entry'].get()) for r in self.resources]
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
                    result = self._solve_simplex(list(product_combo), b_vector, num_resources, opt_type)
                    results_list.append(result)

                if not results_list:
                    messagebox.showinfo("Brak wyników", "Nie udało się obliczyć żadnej kombinacji.")
                    return

                valid_results = [r for r in results_list if r['status'] in ('optimal', 'unbounded')]
                if not valid_results:
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

                result = self._solve_simplex(all_products_data, b_vector, num_resources, opt_type)
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