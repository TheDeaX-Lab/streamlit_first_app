import streamlit as st
import pandas as pd
import graphviz
from functools import reduce
from fractions import Fraction
import numpy as np

def _sum(t1, t2):
    return f"{t1} + {t2}"

def _mul(t1, t2):
    return f"{t1} \cdot {t2}"

nline = "\\\\"

st.title("Строим матрицу определенной размерности")
st.subheader("Настройки")
n = st.number_input("Размерность матрицы", 2, 10)
chk = st.checkbox("Состояние нумеруются от единицы?")
length_digits = st.number_input("Цифр после запятой", min_value=1, max_value=10)
cols = [f"S{i + chk}" for i in range(n)]
st.subheader("Задание значений матрицы перехода состояний")

rawlistvalues = []
for i in range(n):
    st.subheader(f"Значения переходов из {cols[i]}")
    rw = []
    for j in range(n):
        rw.append(
            Fraction(
                st.number_input(f"Вероятность перехода в {cols[j]}. Значение вводите в соответствии с умножением вероятности на {10**length_digits}", value=10**length_digits if i == j else 0, min_value=0, max_value=10**length_digits, key=(i, j)), (10**length_digits))
        )
    rawlistvalues.append(rw)

valss = np.asarray([[rawlistvalues[i][j] for j in range(n)] for i in range(n)])
df = pd.DataFrame(valss, columns=cols, index=cols, dtype="float64")
# editable_tb = st.data_editor(df)
st.subheader("Предпросмотр полученной матрицы")
st.dataframe(df)
dff = df.T.sum()

enable_calculate = st.checkbox("Включить вычисления?")
if enable_calculate:
    if np.all(np.sum(valss, axis=1) == 1):
        st.title("Проделанные вычисления")
        st.header("Граф состояний")
        vals = df.values
        graph = graphviz.Graph()
        for i in range(n):
            for j in range(n):
                val = vals[i][j]
                if val > 0:
                    graph.edge(f"{cols[i]}", f"{cols[j]}", f"{val}")
        st.graphviz_chart(graph)

        pjcols = [f"P_{i + chk}" for i in range(n)]
        pijcols = [["p_{" + f"{i + chk}{j + chk}" + "}" for j in range(n)] for i in range(n)]
        
        st.header("Общая формула в зависимости от шага")
        for j in range(n):
            st.latex(pjcols[j] + "(t) = " + reduce(_sum, [_mul(pjcols[i] + "(t - 1)", pijcols[i][j]) for i in range(n)])
            )
        st.header("Настройки для расчета с нулевого шага")
        nsteps = st.number_input("Количество шагов", min_value=0, max_value=10)
        st.header("Вероятности при шаге 0")
        rw = []
        for j in range(n):
            rw.append(
                Fraction(st.number_input(f"Исходная вероятность для состояния {cols[j]}. Значение вводите в соответствии с умножением вероятности на {10**length_digits}", value=10**length_digits if j == 0 else 0, min_value=0, max_value=10**length_digits,
                                key=(j)), 10**length_digits)
            )
        rw = np.asarray(rw)
        steps = [rw]
        
        for step in range(nsteps):
            steps.append(valss @ steps[-1])
        Pdf = pd.DataFrame(steps, columns=cols, dtype="float64")
        Pdf.index = [f"k={i}" for i in range(len(steps))]
        Pdf.columns.name = "Вероятность Pi(k)"
        Pdf.index.name = "Количество шагов k"
        Pdf["SUM"] = Pdf.apply(sum, axis=1)
        st.header("Итоги вычислений")
        st.dataframe(Pdf)

        st.header("Общая формула стационарного уравнения")
        for j in range(n):
            st.latex(pjcols[j] + " = " + reduce(_sum, [_mul(pjcols[i], pijcols[i][j]) for i in range(n)])
            )

        st.header("СЛАУ для решения стационарного уравнения")
        for j in range(n):
            st.latex(pjcols[j] + " = " + (reduce(_sum, [_mul(pjcols[i], f"{vals[i][j]}") for i in range(n) if vals[i][j] != 0]) if sum(vals[:, j]) != 0 else "0")
            )
        
        st.header("Попытка автоматического решения СЛАУ (Может не сработать)")
        
        ns = valss.copy()
        for i in range(len(ns)):
            ns[i][i] -= 1
        for i in range(len(ns)):
            if ns[i:, i].sum() == 0: continue
            if ns[i, i] == 0:
                for k in range(i + 1, len(ns)):
                    if ns[k][i] != 0:
                        break
                ns[i], ns[k] = ns[k], ns[i]
            ns[i] = ns[i] / ns[i][i]
            for j in range(i + 1, len(ns)):
                ns[j] -= ns[i] * ns[j][i]
        for i in reversed(list(range(1, len(ns) - 1))):
            if ns[i][i] == 0: continue
            ns[i] = ns[i] / ns[i][i]
            for j in range(i):
                ns[j] -= ns[i] * ns[j][i]
        import sympy
        Ps = {
            k: sympy.Symbol(k) for k in pjcols
        }
        result_values = sympy.solve([
            sympy.Eq(
                sympy.Add(
                    *[
                        Ps[pjcols[k]] * v
                        for k, v in enumerate(row)
                    ]
                ), 0
            )
            for row in ns
            if sum(map(abs, row)) != 0
        ] + [sympy.Eq(sympy.Add(*Ps.values()), 1)], Ps)
        if len(result_values) > 0:
            for k in range(len(ns)):
                if Ps[pjcols[k]] in result_values:
                    st.latex(pjcols[k] + " = " + str(result_values[Ps[pjcols[k]]]))
                else:
                    st.latex(pjcols[k] + " = ХЗ")
            
            def tauj(j):
                return (q / (1 - valss[j][j])) if valss[j][j] != 1 else float("inf")

            def taurj(j):
                return tauj(j)* (
                    (2 - valss[j][j] - result_values[Ps[pjcols[j]]])
                    /
                    result_values[Ps[pjcols[j]]]
                ) if Ps[pjcols[k]] in result_values else "ХЗ"

            def tj(j):
                return (tauj(j) + taurj(j)) if Ps[pjcols[k]] in result_values else "ХЗ"

            st.header("Опции для расчета с квантилем времени")
            d = st.number_input("Сколько цифр после запятой", min_value=0, max_value=10, value=0)
            q = st.number_input("Квантиль времени с учетом делителя выше", min_value=1, value=1) / (10**d)

            df = pd.DataFrame(cols, columns=["Состояние"])
            df["t_j"] = [tauj(i) for i in range(n)]
            df["t_rj"] = [taurj(i) for i in range(n)]
            df["T_j"] = [tj(i) for i in range(n)]
            df.set_index("Состояние", inplace=True)
            df.columns.name = "Характеристики в секундах"
            df.rename(columns={"t_j": "Время пребывания внутри состояния", "t_rj": "Время пребывания вне состояния", "T_j": "Общее время"}, inplace=True)

            st.dataframe(df)
        else:
            st.write("Увы решение не удалось")
    else:
        st.write("Таблица состояний составлена некорректно, сумма по строке должна быть 1")
        st.dataframe(pd.DataFrame([*zip(dff, np.sum(valss, axis=1) == 1)], columns=["Сумма по строке", "Допущен?"], index=cols))