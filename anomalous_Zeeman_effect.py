# author:       0.382 (https://github.com/0382)
# enviroment:   win10, python3.7.4
# require:      numpy, matplotlib
# description:  反常塞曼效应作图
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist
from fractions import Fraction
import json

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']

# L=0时,用符号 `S` 表示,...
L_marker = {
    0 : 'S',
    1 : 'P',
    2 : 'D',
    3 : 'F',
    4 : 'G'
}

def LS_g_factor(L:Fraction, S:Fraction, J:Fraction)->Fraction:
    '''LS 耦合的 g 因子, Jj 耦合太麻烦了,不写了'''
    return 1 + (J*(J+1) - L*(L+1) + S*(S+1)) / (2*J*(J+1))

def get_m_list(J:Fraction)->list:
    '''给定 J 值, 返回可能的 M 值列表'''
    N:int = int(J*2+1)
    return [-J + i for i in range(N)]

def anomalous_Zeeman_effect_transition(LSJ1:list, LSJ2:list)->dict:
    '''
    获得反常塞曼效应跃迁结果
    Parameters
    ----------
    LSJ1: 上能级 L1, S1, J1 列表
    LSJ2: 下能级 L2, S2, J2 列表

    Returns
    -------
    返回形如
    result: dict = {
        'pi': [
            (m1, m2, intensity),
            ...
        ]
        'sigma': [
            (m1, m2, intensity),
            ....
        ]
    }
    m1, m2 为上下能级的 m 量子数
    intensity 为子谱线相对强度
    pi 和 sigma 区别 pi 光 (M->M) 和 sigma 光 (M->M \pm 1)
    如果不发生跃迁, 返回空字典
    '''
    LSJ1 = [Fraction(__) for __ in LSJ1]
    LSJ2 = [Fraction(__) for __ in LSJ2]
    L1, S1, J1 = LSJ1
    L2, S2, J2 = LSJ2
    # 选择定则
    assert((J1-J2) in [0, 1, -1])
    m1_list = get_m_list(J1)
    m2_list = get_m_list(J2)
    N1 = int(2*J1+1)
    N2 = int(2*J2+1)
    result = {'pi':[], 'sigma': []}
    # 如果 J -> J
    if J1 == J2:
        if J1 == 0:
            return {}
        # pi 光
        for i in range(N1):
            if m1_list[i] == 0:
                continue
            result['pi'].append((m1_list[i], m2_list[i], m1_list[i]**2))
        # sigma 光
        for i in range(N1):
            m1 = m1_list[i]
            if i == 0:
                result['sigma'].append((m1, m2_list[i+1], Fraction(1,4)*(J1+m1+1)*(J1-m1)))
            elif i == N1-1:
                result['sigma'].append((m1, m2_list[i-1], Fraction(1,4)*(J1-m1+1)*(J1+m1)))
            else:
                result['sigma'].append((m1, m2_list[i+1], Fraction(1,4)*(J1+m1+1)*(J1-m1)))
                result['sigma'].append((m1, m2_list[i-1], Fraction(1,4)*(J1-m1+1)*(J1+m1)))
    # 如果 J -> J+1
    elif J1 == J2-1:
        # pi 光
        for i in range(N1):
            assert(m1_list[i] == m2_list[i+1])
            result['pi'].append((m1_list[i], m2_list[i+1], (J1+1)**2 - m1_list[i]**2))
        # sigma 光
        for i in range(N1):
            m1 = m1_list[i]
            result['sigma'].append((m1, m2_list[i], Fraction(1,4)*(J1-m1+1)*(J1-m1+2)))
            result['sigma'].append((m1, m2_list[i+2], Fraction(1,4)*(J1+m1+1)*(J1+m1+2)))
    # 如果 J -> J-1
    elif J1 == J2+1:
        # pi 光
        for i in range(N2):
            assert(m2_list[i] == m1_list[i+1])
            result['pi'].append((m1_list[i+1], m2_list[i], J1**2 - m1_list[i+1]**2))
        # sigma 光
        for i in range(N2):
            m2 = m2_list[i]
            result['sigma'].append((m1_list[i], m2, Fraction(1,4)*(J1-m2+1)*(J1-m2)))
            result['sigma'].append((m1_list[i+2], m2, Fraction(1,4)*(J1+m2+1)*(J1+m2)))
    else:
        result = {}
    if result != {}:
        result['pi'].sort()
        result['sigma'].sort()
    return result

def azet_result_to_str(result:dict)->dict:
    '''因为结果含有分数,不容易阅读,转化为字符格式便于阅读,
    分数的结果后面还有用,因此写一个单独的转换函数'''
    sr = {'pi': [], 'sigma': []}
    for v in result['pi']:
        sr['pi'].append(f'{v[0]}, {v[1]}, {v[2]}')
    for v in result['sigma']:
        sr['sigma'].append(f'{v[0]}, {v[1]}, {v[2]}')
    return json.dumps(sr, indent=2)

def anomalous_Zeeman_effect_split(LSJ1:list, LSJ2:list, filename:str=None):
    '''对反常塞曼效应做跃迁能级图
    Parameters:
    LSJ1: 上能级的 l, s, j 值，通过列表传入，如 [l, 1/2, 3/2]
    LSJ2: 下能级的 l, s, j 值
    filename: 保存为文件,如果参数没有指定,不保存为文件,会调用`plot.show()`展示出来

    只是对 LS 耦合做了一些计算, Jj耦合太过于复杂, 就不写了
    '''
    LSJ1 = [Fraction(__) for __ in LSJ1]
    LSJ2 = [Fraction(__) for __ in LSJ2]
    L1, S1, J1 = LSJ1
    L2, S2, J2 = LSJ2
    g1 = LS_g_factor(*LSJ1)
    g2 = LS_g_factor(*LSJ2)
    print(f'g1 = {g1}, g2 = {g2}')

    plt.figure()
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.grid()
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    level_1 = 25        #作图的一些定标常数
    level_2 = 75
    vertical_1 = 0
    vertical_2 = 30
    vertical_3 = 40
    vertical_4 = 75
    vertical_5 = 76
    vertical_6 = 85
    swell_factor = 5    # 能级图的相对间距，如果图不好看，试着调节一下这个参数

    m12mg = lambda m1 : float(m1 * g1 * swell_factor)
    m22mg = lambda m2 : float(m2 * g2 * swell_factor)

    # 能级符号
    plt.text(vertical_1, level_2, r"$^{%s}\mathrm{%s}_{%s}$" % (f'{2*S1+1}', L_marker[int(L1)], f'{J1}'))
    plt.text(vertical_1, level_1, r"$^{%s}\mathrm{%s}_{%s}$" % (f'{2*S2+1}', L_marker[int(L2)], f'{J2}'))
    # 原能级线
    plt.plot([vertical_1+8, vertical_2], [level_1, level_1], '-', color='black', linewidth=1)
    plt.plot([vertical_1+8, vertical_2], [level_2, level_2], '-', color='black', linewidth=1)

    M1_list = get_m_list(J1)
    M2_list = get_m_list(J2)
    lv_1_list = [m12mg(__) + level_2 for __ in M1_list]
    lv_2_list = [m22mg(__) + level_1 for __ in M2_list]
    # 分裂能级线
    for lv in lv_1_list:
        # 原能级到分裂能级的虚线
        plt.plot([vertical_2, vertical_3], [level_2, lv], linestyle=(8,(5,5)), color='black', linewidth=1)
        # 分裂能级线
        plt.plot([vertical_3, vertical_4], [lv, lv], '-', color='black', linewidth=1)
    for lv in lv_2_list:
        plt.plot([vertical_2, vertical_3], [level_1, lv], linestyle=(8,(5,5)), color='black', linewidth=1)
        plt.plot([vertical_3, vertical_4], [lv, lv], '-', color='black', linewidth=1)

    # M 值
    plt.text(vertical_5, lv_1_list[-1]+5, '$M$')
    plt.text(vertical_6, lv_1_list[-1]+5, '$Mg$')
    for i,m in enumerate(M1_list):
        plt.text(vertical_5, lv_1_list[i], f'{m}')
        plt.text(vertical_6, lv_1_list[i], f'{m*g1}')
    for i,m in enumerate(M2_list):
        plt.text(vertical_5, lv_2_list[i], f'{m}')
        plt.text(vertical_6, lv_2_list[i], f'{m*g2}')
    
    # 箭头
    r = anomalous_Zeeman_effect_transition(LSJ1, LSJ2)
    dx = (vertical_4 - vertical_3) / (len(r['pi']) + len(r['sigma']) + 1)
    arr_x = vertical_3 + dx/2
    for i,spec in enumerate(r['pi']):
        y1 = m12mg(spec[0]) + level_2
        y2 = m22mg(spec[1]) + level_1
        plt.annotate("", xy=[arr_x, y2-1], xytext=[arr_x, y1], arrowprops=dict(arrowstyle="->"))
        plt.text(arr_x, y1+1, r"$\pi_%d$"%(i+1))
        arr_x += dx
    for i,spec in enumerate(r['sigma']):
        y1 = m12mg(spec[0]) + level_2
        y2 = m22mg(spec[1]) + level_1
        plt.annotate("", xy=[arr_x, y2-1], xytext=[arr_x, y1], arrowprops=dict(arrowstyle="->"))
        plt.text(arr_x, y1+1, r"$\sigma_%d$"%(i+1))
        arr_x += dx

    # 有无磁场的标记
    plt.text(vertical_1+15, level_2+3, '无磁场')
    plt.text(vertical_3+15, lv_1_list[-1]+3, '有磁场')

    if(filename == None):
        plt.show()
    else:
        plt.savefig(filename)

def anomalous_Zeeman_effect_intensity(LSJ1:list, LSJ2:list, filename:str=None):
    '''对反常塞曼效应作相对强度分布图
    Parameters:
    LSJ1: 上能级的 l, s, j 值，通过列表传入，如 [l, 1/2, 3/2]
    LSJ2: 下能级的 l, s, j 值
    filename: 保存为文件,如果参数没有指定,不保存为文件,会调用`plot.show()`展示出来
    '''
    LSJ1 = [Fraction(__) for __ in LSJ1]
    LSJ2 = [Fraction(__) for __ in LSJ2]
    azet = anomalous_Zeeman_effect_transition(LSJ1, LSJ2)
    g1 = LS_g_factor(*LSJ1)
    g2 = LS_g_factor(*LSJ2)
    pi = [[g1*m1 - g2*m2, intensity] for m1,m2,intensity in azet['pi']]
    sigma = [[g1*m1 - g2*m2, intensity] for m1,m2,intensity in azet['sigma']]
    positions = [pos for pos, intensity in pi]+[pos for pos, intensity in sigma]
    
    fig = plt.figure()
    ax = axisartist.Subplot(fig, 1,1,1)
    fig.add_axes(ax)
    max_pos = float(max(positions))
    max_intensity = float(max([intensity for pos, intensity in pi]+[pos for pos, intensity in sigma]))
    plt.xlim(-max_pos - 0.2, max_pos + 0.2)
    plt.ylim(-max_intensity, max_intensity)
    ax.axis[:].set_visible(False)
    ax.axis['x'] = ax.new_floating_axis(0,0)
    ax.axis['x'].set_axisline_style("->", size = 2.0)
    ax.axis["x"].set_axis_direction('top')
    ax.set_xticks([float(pos) for pos in positions])
    ax.set_xticklabels([str(pos) for pos in positions])
    vertical = lambda x, y : plt.plot([x, x], [0, float(y)], '-', color='black', linewidth=1)
    text = lambda x, y : plt.text(x, float(y), f'{y}')
    for i,spectrum in enumerate(pi):
        vertical(spectrum[0], spectrum[1])
        text(spectrum[0], spectrum[1])
        plt.text(spectrum[0], spectrum[1]/2, r'$\pi_%d$'%(i+1))
    for i,spectrum in enumerate(sigma):
        vertical(spectrum[0], -spectrum[1])
        text(spectrum[0], -spectrum[1])
        plt.text(spectrum[0], -spectrum[1]/2, r'$\sigma_%d$'%(i+1))
    plt.text(max_pos, -max_intensity/2, r'$\sigma$光')
    plt.text(max_pos, max_intensity/2, r'$\pi$光')
    if filename == None:
        plt.show()
    else:
        plt.savefig(filename)

if __name__ == "__main__":
    # 例如
    # 钠灯 589.0 nm, 上能级 2P3/2, 下能级 2S1/2
    LSJ1 = [1, 1/2, 3/2]
    LSJ2 = [0, 1/2, 1/2]
    # 汞灯 546.1 nm, 上能级 3S1, 下能级 3P2
    LSJ1 = [0, 1, 1]
    LSJ2 = [1, 1, 2]
    anomalous_Zeeman_effect_split(LSJ1, LSJ2, "images/Hg-split.png")
    r:dict = anomalous_Zeeman_effect_transition(LSJ1, LSJ2)
    print(azet_result_to_str(r))
    anomalous_Zeeman_effect_intensity(LSJ1, LSJ2, "images/Hg-intensity.png")