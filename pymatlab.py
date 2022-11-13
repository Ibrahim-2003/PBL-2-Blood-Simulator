
text = 'Elv_s Elv_d Cpas Cpat Cpvn Erv_s Erv_d Csas Csat Csvn Rpas Lpas Rpat Rpar Rpcp Rpvn Lsas'
text+= ' Rsas Rsat Lsat Rsar Rscp Rsvn CQao CQmi Plv_0'
text+= ' Vlv_0 Ela_max Ela_min Pla_0 Vla_0 CQpo CQti'
text+= ' Erv_s Erv_d Prv_0 Vrv_0 Era_max Era_min Pra_0'
text+=' Vra_0 Kst_la Kst_lv Kf_sav Ke_sav Msav Asav'
text+=" Kst_ra Kst_rv Kf_pav Ke_pav Mpav Apav Kp_ao"
text+= ' Kf_ao Kp_mi Kf_mi Kp_sv Kf_sv Kp_po Kp_ti'
text+= ' Kf_ti Kp_pv Kf_pv theta_ao_max theta_po_max'
text+=' theta_ti_max theta_mi_max Rsm Rbrain Rkid'
text+=' Rliv Rsp Rint T Tpwb Tpww Ts1 Ts2 Lpat Kf_po'

texts = [i for i in text.split(' ')]
test = []
for i in range(len(texts)):
    if texts[i] in test:
        print(texts[i])
    test.append(texts[i])

