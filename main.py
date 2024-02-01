import streamlit as st
import pandas as pd
import graphviz
from functools import reduce

st.title("Строим матрицу определенной размерности")

n = st.number_input("Размерность матрицы", 2, 10)
chk = st.checkbox("Состояние нумеруются от единицы?")

cols = [f"S{i + chk}" for i in range(n)]
divider = (n * (n - 1)) // 2
df = pd.DataFrame([[i / divider for i in range(n)] for _ in range(n)], columns=cols, index=cols, dtype="float64")
editable_tb = st.data_editor(df)
dff = editable_tb.T.sum()

enable_calculate = st.checkbox("Включить вычисления?")
if enable_calculate:
    if (dff == 1).all():
        st.title("Граф состояний")
        graph = graphviz.Graph()
        vals = editable_tb.values
        for i in range(n):
            for j in range(n):
                val = vals[i][j]
                if val > 0:
                    graph.edge(f"{cols[i]}", f"{cols[j]}", f"{val}")
        st.graphviz_chart(graph)

        pjcols = [f"P_{i + chk}" for i in range(n)]
        pijcols = [["p_{" + f"{i + chk}{j + chk}" + "}" for j in range(n)] for i in range(n)]

        
        def _sum(t1, t2):
            return f"{t1} + {t2}"

        def _mul(t1, t2):
            return f"{t1} \cdot {t2}"
        
        nline = "\\\\"
        
        st.title("Общая формула в зависимости от шага")
        for j in range(n):
            st.latex(pjcols[j] + "(t) = " + reduce(_sum, [_mul(pjcols[i] + "(t - 1)", pijcols[i][j]) for i in range(n)])
            )

        st.title("Общая формула стационарного")
        for j in range(n):
            st.latex(pjcols[j] + " = " + reduce(_sum, [_mul(pjcols[i], pijcols[i][j]) for i in range(n)])
            )

        st.title("СЛАУ для решения стационарного уравнения")
        for j in range(n):
            st.latex(pjcols[j] + " = " + (reduce(_sum, [_mul(pjcols[i], f"{vals[i][j]}") for i in range(n) if vals[i][j] != 0]) if sum(vals[:, j]) != 0 else "0")
            )
    else:
        st.write("Таблица состояний составлена некорректно, сумма по строке должна быть 1")
        st.dataframe(dff)