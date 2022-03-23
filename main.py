from iapws import IAPWS97 as WSP
import math as M
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color='r'):
    if not isinstance(wsp_point, list):
        iso_bar_0_s = np.arange(wsp_point.s + min_s, wsp_point.s + max_s, step_s).tolist()
        iso_bar_0_h = [WSP(P=wsp_point.P, s=i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s + min_s, wsp_point[1].s + max_s, step_s).tolist()
        iso_bar_0_h = [WSP(P=wsp_point[1].P, s=i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s, iso_bar_0_h, color)
d = 1.1  # m
p_0 = 26.5  # МПа
T_0 = 873.15  # K
n = 60  # Гц
G_0 = 425.7379  # кг/с
H_0 = 105  # кДж/кг
rho = 0.05
l_1 = 0.0265  # м
alpha_1 = 12  # град
b_1 = 0.06  # м
Delta = 0.003  # м
b_2 = 0.03  # м
kappa_vs = 0  # коэф исп вых скорости
st.header("Курсовая работа вторая часть. Юричковский.К.А ТФэ-01-19 Вариант 23")
st.write("https://github.com/KirillYur/Turbine-2")
def callculate_optimum(d, p_0, T_0, n, G_0, H_0, rho, l_1, alpha_1, b_1, Delta, b_2, kappa_vs):
    u = M.pi * d * n
    point_0 = WSP(P=p_0, T=T_0)
    H_0s = H_0 * (1 - rho)
    H_0r = H_0 * rho
    h_1t = point_0.h - H_0s
    point_1t = WSP(h=h_1t, s=point_0.s)
    c_1t = (2000 * H_0s) ** 0.5
    M_1t = c_1t / point_1t.w
    mu_1 = 0.982 - 0.005 * (b_1 / l_1)
    F_1 = G_0 * point_1t.v / mu_1 / c_1t
    el_1 = F_1 / M.pi / d / M.sin(M.radians(alpha_1))
    e_opt = 5 * el_1 ** 0.5
    if e_opt > 0.85:
        e_opt = 0.85
    l_1 = el_1 / e_opt

    fi_1 = 0.98 - 0.008 * (b_1 / l_1)
    c_1 = c_1t * fi_1
    alpha_1 = M.degrees(M.asin(mu_1 / fi_1 * M.sin(M.radians(alpha_1))))
    w_1 = (c_1 ** 2 + u ** 2 - 2 * c_1 * u * M.cos(M.radians(alpha_1))) ** 0.5
    betta_1 = M.degrees(M.atan(M.sin(M.radians(alpha_1)) / (M.cos(M.radians(alpha_1)) - u / c_1)))
    Delta_Hs = c_1t ** 2 / 2 * (1 - fi_1 ** 2)
    h_1 = h_1t + Delta_Hs * 1e-3
    point_1 = WSP(P=point_1t.P, h=h_1)
    h_2t = h_1 - H_0r
    point_2t = WSP(h=h_2t, s=point_1.s)
    w_2t = (2 * H_0r * 1e3 + w_1 ** 2) ** 0.5
    l_2 = l_1 + Delta
    mu_2 = 0.965 - 0.01 * (b_2 / l_2)
    M_2t = w_2t / point_2t.w
    F_2 = G_0 * point_2t.v / mu_2 / w_2t
    betta_2 = M.degrees(M.asin(F_2 / (e_opt * M.pi * d * l_2)))
    point_1w = WSP(h=point_1.h + w_1 ** 2 / 2 * 1e-3, s=point_1.s)

    psi = 0.96 - 0.014 * (b_2 / l_2)
    w_2 = psi * w_2t
    c_2 = (w_2 ** 2 + u ** 2 - 2 * u * w_2 * M.cos(M.radians(betta_2))) ** 0.5
    alpha_2 = M.degrees(M.atan(M.sin(M.radians(betta_2)) / (M.cos(M.radians(betta_2)) - u / w_2)))
    if alpha_2 < 0:
        alpha_2 = 180 + alpha_2
    Delta_Hr = w_2t ** 2 / 2 * (1 - psi ** 2)
    h_2 = h_2t + Delta_Hr * 1e-3
    point_2 = WSP(P=point_2t.P, h=h_2)
    Delta_Hvs = c_2 ** 2 / 2
    E_0 = H_0 - kappa_vs * Delta_Hvs
    etta_ol1 = (E_0 * 1e3 - Delta_Hs - Delta_Hr - (1 - kappa_vs) * Delta_Hvs) / (E_0 * 1e3)
    etta_ol2 = (u * (c_1 * M.cos(M.radians(alpha_1)) + c_2 * M.cos(M.radians(alpha_2)))) / (E_0 * 1e3)
    return etta_ol2, alpha_2

d = [i * 1e-2 for i in list(range(90, 111, 1))]
alpha1 = []
eta = []
ucf = []
fighs = plt.figure()
for i in d:
    ucf_1 = M.pi * i * n / (2000 * H_0) ** 0.5
    ucf.append(ucf_1)

    eta_ol, alpha = callculate_optimum(i, p_0, T_0, n, G_0, H_0, rho, l_1, alpha_1, b_1, Delta, b_2, kappa_vs)
    alpha1.append(alpha)
    eta.append(eta_ol)
plt.plot(ucf, eta)
plt.ylabel('eta')
plt.xlabel('u_cf')
st.pyplot(fighs)
def frange(x, y, jump):
  while x < y:
    yield x
    x += jump
df = pd.DataFrame({
    "d, м" : list(frange(0.9, 1.11, 0.01)),
    "eta_ol" : (eta),
    "alpha" : (alpha1),
    "U_cf" : (ucf)})
df
d = 1.1
u = M.pi * d * n
with st.expander("See explanation"):
    st.write(f'u = {u:.2f} м/с')
    point_0 = WSP(P=p_0, T=T_0)
    st.write(f'h_0 = {point_0.h:.2f} кДж/кг')
    st.write(f's_0 = {point_0.s:.4f} кДж/(кг*К)')
    H_0s = H_0 * (1 - rho)
    H_0r = H_0 * rho
    h_1t = point_0.h - H_0s
    st.write(f'h_1т = {h_1t:.2f} кДж/кг')
    point_1t = WSP(h=h_1t, s=point_0.s)
    c_1t = (2000 * H_0s) ** 0.5
    st.write(f'c_1т = {c_1t:.2f} м/с')
    M_1t = c_1t / point_1t.w
    st.write(f'M_1т = {M_1t:.2f}')
    mu_1 = 0.982 - 0.005 * (b_1 / l_1)
    F_1 = G_0 * point_1t.v / mu_1 / c_1t
    st.write(f'F_1 = {F_1:.4f} м^2')
    el_1 = F_1 / M.pi / d / M.sin(M.radians(alpha_1))
    st.write(f'el_1 = {el_1:.4f} м')
    e_opt = 6 * el_1 ** 0.5
if e_opt > 0.85:
    e_opt = 0.85
    l_1 = el_1 / e_opt
    st.write(f'l_1 = {l_1:.4f} м')

def plot_hs_nozzle_t(x_lim, y_lim):
    plt.plot([point_0.s, point_1t.s], [point_0.h, point_1t.h], 'ro-')
    iso_bar(point_0, -0.02, 0.02, 0.001, 'c')
    iso_bar(point_1t, -0.02, 0.02, 0.001, 'y')
    plt.xlim(x_lim)
    plt.ylim(y_lim)
plot_hs_nozzle_t([6.1, 6.5], [3300, 3600])
plt.ylabel('H кДж/кг')
plt.xlabel('S кДж/кг*К')

if alpha_1 <= 10:
    NozzleBlade = 'C-90-09A'
    t1_ = 0.78
    b1_mod = 6.06
    f1_mod = 3.45
    W1_mod = 0.471
    alpha_inst1 = alpha_1 - 12.5 * (t1_ - 0.75) + 20.2
elif 10 < alpha_1 <= 13:
    NozzleBlade = 'C-90-12A'
    t1_ = 0.78
    b1_mod = 5.25
    f1_mod = 4.09
    W1_mod = 0.575
    alpha_inst1 = alpha_1 - 10 * (t1_ - 0.75) + 21.2
elif 13 < alpha_1 <= 16:
    NozzleBlade = 'C-90-15A'
    t1_ = 0.78
    b1_mod = 5.15
    f1_mod = 3.3
    W1_mod = 0.45
    alpha_inst1 = alpha_1 - 16 * (t1_ - 0.75) + 23.1
else:
    NozzleBlade = 'C-90-18A'
    t1_ = 0.75
    b1_mod = 4.71
    f1_mod = 2.72
    W1_mod = 0.333
    alpha_inst1 = alpha_1 - 17.7 * (t1_ - 0.75) + 24.2

#st.write('Тип профиля:', NozzleBlade)
#st.write(f'Оптимальный относительный шаг t1_ = {t1_}')
z1 = (M.pi * d) / (b_1 * t1_)
z1 = int(z1)
#if z1 % 2 == 0:

#    st.write(f'z1 = {z1}')
#else:
#    z1 = z1 + 1

#    st.write(f'z1 = {z1}')
t1_ = (M.pi * d) / (b_1 * z1)
Ksi_1_ = (0.021042 * b_1 / l_1 + 0.023345) * 100
k_11 = 7.18977510 * M_1t ** 5 - 26.94497258 * M_1t ** 4 + 39.35681781 * M_1t ** 3 - 26.09044664 * M_1t ** 2 + 6.75424811 * M_1t + 0.69896998
k_12 = 0.00014166 * 90 ** 2 - 0.03022881 * 90 + 2.61549380
k_13 = 13.25474043 * t1_ ** 2 - 20.75439502 * t1_ + 9.12762245
Ksi_1 = Ksi_1_ * k_11 * k_12 * k_13

fi_1 = M.sqrt(1 - Ksi_1 / 100)

#st.write(f'mu_1 = {mu_1}')
#st.write(f'fi_1 = {fi_1}')

alpha_1 = 12
c_1 = c_1t * fi_1

alpha_1 = M.degrees(M.asin(mu_1 / fi_1 * M.sin(M.radians(alpha_1))))

w_1 = (c_1 ** 2 + u ** 2 - 2 * c_1 * u * M.cos(M.radians(alpha_1))) ** 0.5

#st.write(f'c_1 = {c_1:.2f} м/с')
#st.write(f'alpha_1 = {alpha_1:.2f} град.')
#st.write(f'w_1 = {w_1}')
c_1u = c_1 * M.cos(M.radians(alpha_1))
c_1a = c_1 * M.sin(M.radians(alpha_1))
w_1u = c_1u - u

#st.write(c_1u, w_1u)
w_1_tr = [0, 0, -w_1u, -c_1a]
c_1_tr = [0, 0, -c_1u, -c_1a]
u_1_tr = [-w_1u, -c_1a, -u, 0]

ax = plt.axes()
ax.arrow(*c_1_tr, head_width=5, length_includes_head=True, head_length=20, fc='r', ec='r')
ax.arrow(*w_1_tr, head_width=5, length_includes_head=True, head_length=20, fc='b', ec='b')
ax.arrow(*u_1_tr, head_width=5, length_includes_head=True, head_length=20, fc='g', ec='g')
plt.text(-2 * c_1u / 3, -3 * c_1a / 4, '$c_1$', fontsize=20)
plt.text(-2 * w_1u / 3, -3 * c_1a / 4, '$w_1$', fontsize=20)
betta_1 = M.degrees(M.atan(M.sin(M.radians(alpha_1)) / (M.cos(M.radians(alpha_1)) - u / c_1)))

Delta_Hs = c_1t ** 2 / 2 * (1 - fi_1 ** 2)


h_1 = h_1t + Delta_Hs * 1e-3

point_1 = WSP(P=point_1t.P, h=h_1)
h_2t = h_1 - H_0r

point_2t = WSP(h=h_2t, s=point_1.s)
w_2t = (2 * H_0r * 1e3 + w_1 ** 2) ** 0.5

l_2 = l_1 + Delta
mu_2 = 0.965 - 0.01 * (b_2 / l_2)

M_2t = w_2t / point_2t.w

F_2 = G_0 * point_2t.v / mu_2 / w_2t

betta_2 = M.degrees(M.asin(F_2 / (e_opt * M.pi * d * l_2)))

point_1w = WSP(h=point_1.h + w_1 ** 2 / 2 * 1e-3, s=point_1.s)
#st.write(f'betta_1 = {betta_1:.2f}')
#st.write(f'Delta_Hs = {Delta_Hs:.2f} Дж/кг')
#st.write(f'h_1 = {h_1:.2f} кДж/кг')
#st.write(f'h_2t = {h_2t:.2f} кДж/кг')
#st.write(f'w_2t = {w_2t:.2f} м/с')
#st.write(f'mu_2 = {mu_2:.2f}')
#st.write(f'M_2t = {M_2t:.2f}')
#st.write(f'F_2 = {F_2:.2f}')
#st.write(f'betta_2 = {betta_2:.2f}')
fig3 = plt.figure()
def plot_hs_stage_t(x_lim, y_lim):
    plot_hs_nozzle_t(x_lim, y_lim)
    plt.plot([point_0.s, point_1.s], [point_0.h, point_1.h], 'bo-')
    plt.plot([point_1.s, point_2t.s], [point_1.h, point_2t.h], 'ro-')
    plt.plot([point_1.s, point_1.s], [point_1w.h, point_1.h], 'ro-')
    iso_bar(point_2t, -0.02, 0.02, 0.001, 'y')
    iso_bar(point_1w, -0.005, 0.005, 0.001, 'c')
plt.ylabel('H кДж/кг')
plt.xlabel('S кДж/кг*К')

plot_hs_stage_t([6.28, 6.35], [3350, 3500])

st.pyplot(fig3)

if betta_2 <= 15:
    RotorBlade = 'P-23-14A'
    t2_ = 0.63
    b2_mod = 2.59
    f2_mod = 2.44
    W2_mod = 0.39
    beta_inst2 = betta_2 - 12.5 * (t2_ - 0.75) + 20.2

elif 15 < betta_2 <= 19:
    RotorBlade = 'P-26-17A'
    t2_ = 0.65
    b2_mod = 2.57
    f2_mod = 2.07
    W2_mod = 0.225
    beta_inst2 = betta_2 - 19.3 * (t2_ - 0.6) + 60


elif 19 < betta_2 <= 23:
    RotorBlade = 'P-30-21A'
    t2_ = 0.63
    b2_mod = 2.56
    f2_mod = 1.85
    W2_mod = 0.234
    beta_inst2 = betta_2 - 12.8 * (t2_ - 0.65) + 58


elif 23 < betta_2 <= 27:
    RotorBlade = 'P-35-25A'
    t2_ = 0.6
    b2_mod = 2.54
    f2_mod = 1.62
    W2_mod = 0.168
    beta_inst2 = betta_2 - 16.6 * (t2_ - 0.65) + 54.3

elif 27 < betta_2 <= 31:
    RotorBlade = 'P-46-29A'
    t2_ = 0.51
    b2_mod = 2.56
    f2_mod = 1.22
    W2_mod = 0.112
    beta_inst2 = betta_2 - 50.5 * (t2_ - 0.6) + 47.1


else:
    RotorBlade = 'P-50-33A'
    t2_ = 0.49
    b2_mod = 2.56
    f2_mod = 1.02
    W2_mod = 0.079
    beta_inst2 = betta_2 - 20.8 * (t2_ - 0.6) + 43.7



z2 = int((M.pi * d) / (b_2 * t2_))

t2_ = (M.pi * d) / (b_2 * z2)
Ksi_2_ = 4.364 * b_2 / l_2 + 4.22
k_21 = -13.79438991 * M_2t ** 4 + 36.69102267 * M_2t ** 3 - 32.78234341 * M_2t ** 2 + 10.61998662 * M_2t + 0.28528786
k_22 = 0.00331504 * betta_1 ** 2 - 0.21323910 * betta_1 + 4.43127194
k_23 = 60.72813684 * t2_ ** 2 - 76.38053189 * t2_ + 24.97876023
Ksi_2 = Ksi_2_ * k_21 * k_22 * k_23

psi = M.sqrt(1 - Ksi_2 / 100)


psi = 0.93

w_2 = psi * w_2t

c_2 = (w_2 ** 2 + u ** 2 - 2 * u * w_2 * M.cos(M.radians(betta_2))) ** 0.5

alpha_2 = M.degrees(M.atan(M.sin(M.radians(betta_2)) / (M.cos(M.radians(betta_2)) - u / w_2)))

Delta_Hr = w_2t ** 2 / 2 * (1 - psi ** 2)

h_2 = h_2t + Delta_Hr * 1e-3
point_2 = WSP(P=point_2t.P, h=h_2)
Delta_Hvs = c_2 ** 2 / 2

E_0 = H_0 - kappa_vs * Delta_Hvs
etta_ol1 = (E_0 * 1e3 - Delta_Hs - Delta_Hr - (1 - kappa_vs) * Delta_Hvs) / (E_0 * 1e3)

etta_ol2 = (u * (c_1 * M.cos(M.radians(alpha_1)) + c_2 * M.cos(M.radians(alpha_2)))) / (E_0 * 1e3)


c_1u = c_1 * M.cos(M.radians(alpha_1))
c_1a = c_1 * M.sin(M.radians(alpha_1))
w_1u = c_1u - u
w_2a = w_2 * M.sin(M.radians(betta_2))
w_2u = w_2 * M.cos(M.radians(betta_2))
c_2u = w_2u + u

w_1_tr = [0, 0, -w_1u, -c_1a]
c_1_tr = [0, 0, -c_1u, -c_1a]
u_1_tr = [-w_1u, -c_1a, -u, 0]
#st.write('Тип профиля:', RotorBlade)
#st.write(f'Оптимальный относительный шаг t2_ = {t2_}')
#st.write(f'z2 = {z2}')
#st.write(f'psi = {psi:.2f}')
#st.write(f'w_2 = {w_2:.2f} м/с')
#st.write(f'c_2 = {c_2:.2f} м/с')
#st.write(f'alpha_2 = {alpha_2:.2f}')
#st.write(f'Delta_Hr = {Delta_Hr:.2f} Дж/кг')
#st.write(f'Delta_Hvs = {Delta_Hvs:.2f} Дж/кг')
#st.write(f'1. etta_ol = {etta_ol1}')
#st.write(f'2. etta_ol = {etta_ol2}')
w_2_tr = [0, 0, w_2u, -w_2a]
c_2_tr = [0, 0, c_2u, -w_2a]
u_2_tr = [c_2u, -w_2a, -u, 0]
fig4 = plt.figure()
ax = plt.axes()
ax.arrow(*c_1_tr, head_width=5, length_includes_head=True, head_length=20, fc='r', ec='r')
ax.arrow(*w_1_tr, head_width=5, length_includes_head=True, head_length=20, fc='b', ec='b')
ax.arrow(*u_1_tr, head_width=5, length_includes_head=True, head_length=20, fc='g', ec='g')
ax.arrow(*c_2_tr, head_width=5, length_includes_head=True, head_length=20, fc='r', ec='r')
ax.arrow(*w_2_tr, head_width=5, length_includes_head=True, head_length=20, fc='b', ec='b')
ax.arrow(*u_2_tr, head_width=5, length_includes_head=True, head_length=20, fc='g', ec='g')
plt.text(-2 * c_1u / 3, -3 * c_1a / 4, '$c_1$', fontsize=20)
plt.text(-2 * w_1u / 3, -3 * c_1a / 4, '$w_1$', fontsize=20)
plt.text(2 * c_2u / 3, -3 * w_2a / 4, '$c_2$', fontsize=20)
plt.text(2 * w_2u / 3, -3 * w_2a / 4, '$w_2$', fontsize=20)
st.pyplot(fig4)
delta_a = 0.0025
z_per_up = 2
mu_a = 0.5
mu_r = 0.75
d_per = d + l_1
delta_r = d_per * 0.001
delta_ekv = 1 / M.sqrt(1 / (mu_a * delta_a) ** 2 + z_per_up / (mu_r * delta_r) ** 2)

xi_u_b = M.pi * d_per * delta_ekv * etta_ol1 / F_1 * M.sqrt(rho + 1.8 * l_2 / d)

Delta_Hub = xi_u_b * E_0


k_tr = 0.0007
Kappa_VS = 0
u = M.pi * d * n
c_f = M.sqrt(2000 * H_0)
ucf = u / c_f
xi_tr = k_tr * d ** 2 / F_1 * ucf ** 3

Delta_Htr = xi_tr * E_0


k_v = 0.065
m = 1
xi_v = k_v / M.sin(M.radians(alpha_1)) * (1 - e_opt) / e_opt * ucf ** 3 * m

i_p = 4
B_2 = b_2 * M.sin(M.radians(beta_inst2))
xi_segm = 0.25 * B_2 * l_2 / F_1 * ucf * etta_ol1 * i_p

xi_parc = xi_v + xi_segm
Delta_H_parc = E_0 * xi_parc

H_i = E_0 - Delta_Hr * 1e-3 - Delta_Hs * 1e-3 - (1 - Kappa_VS) * Delta_Hvs * 1e-3 - Delta_Hub - Delta_Htr - Delta_H_parc

eta_oi = H_i / E_0

N_i = G_0 * H_i

#st.write("""Эквивалентный зазор в уплотнении по бандажу (периферийном)
#            delta_ekv = %.3f мм""" % (delta_ekv * 1000))
#st.write("""Относительные потери от утечек через бандажные уплотнения
#                    xi_u_b = %.3f""" % xi_u_b)
#st.write("""Абсолютные потери от утечек через периферийное уплотнение ступени
#                     Delta_Hub = %.3f кДж/кг""" % Delta_Hub)
#st.write("""Определяем u/c_ф для ступени
#                     U/c_ф = %.3f""" % ucf)
#st.write("""Относительные потери от трения диска
#                     xi_tr = %.5f""" % xi_tr)
#st.write("""Абсолютные потери от трения диска
#                     Delta_Htr = %.3f кДж/кг""" % Delta_Htr)
#st.write("""Относительные вентиляционные потери""", xi_v)
#st.write("""Относительные сегментные потери""", xi_segm)
#st.write("""Использованный теплоперепад ступени
#                    H_i = %.3f кДж/кг""" % H_i)
st.subheader("""Внутренний относительный КПД ступени  
                 η_oi = %.3f """ % eta_oi)
st.subheader("""Внутреняя мощность ступени  
                     N_i = = %.2f кВт""" % N_i)

