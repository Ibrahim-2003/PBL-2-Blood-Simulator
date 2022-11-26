% Copy Catting the Paper

function reverse_engineering

% Known Variables
syms Csa Csp Cep Csv Cev Cpa Cpp Cpv Vusa Vusp Vuep Vusv Vuev Vupa Vupp Vupv...
    Rsa Rsp Rep Rsv Rev Rpa Rpp Rpv Lsa Lpa Cla Vula Rla P0lv kElv Vulv...
    Emaxlv kRlv Cra Vura Rra P0rv kErv Vurv Emaxrv kRrv Pn ka fes_8 fev_0...
    GEmaxlv GEmaxrv GRsp GRep GVusv GVuev GTs GTv fmin tauz fes0 fesmin...
    fev_8 fcs0 tauEmaxlv tauRsp tauRep tauVusv tauVuev tauTs tauTv DEmaxlv...
    DEmaxrv DRsp DRep DVusv DVuev DTs DTv fmax taup kes kev Emaxlv0...
    Emaxrv0 Rsp0 Rep0 Vusv0 Vuev0 T0 I(t) Vt0 Tsys0 ksys0 t0

% Unknown Variables
syms Ppa(t) Fpa(t) Ppp(t) Ppv(t) Psa(t) For(t) Fol(t) Fsa(t) Psp(t) Pev(t) Psv(t) Pra(t) Pla(t) Vlv(t) Vu(t) Vrv(t) Fil(t)...
    Fir(t) Rlv(t) Plv(t) Pmaxlv(t) phi(t) Tsys(t) T(t) epsilon(t) u(t) Rrv(t) Prv(t) Pmaxrv(t) Ptil(t) Pcs(t)...
    fcs(t) fes(t) fev(t) sigma_Emaxrv(t) sigma_Emaxlv(t) sigma_Rsp(t) sigma_Rep(t) sigma_Vusv(t)...
    sigma_Vuev(t) delta_Emaxrv(t) delta_Emaxlv(t) delta_Rsp(t) delta_Rep(t) delta_Vusv(t)...
    delta_Vuev(t) sigma_Ts(t) delta_Ts(t) sigma_Tv(t) delta_Tv(t) Vt(t) tau Gb

eqn1 = diff(Ppa(t),t) == (For(t) - Fpa(t))/Cpa;
eqn2 = diff(Fpa(t),t) == (Ppa(t) - Ppp(t) - Rpa*Fpa(t))/Lpa;
eqn3 = diff(Psa(t),t) == (Fol(t) - Fsa(t))/Csa;
eqn4 = diff(Ppp(t),t) == (Fpa(t) - ((Ppp(t) - Ppv(t))/Rpp))/Cpp;
eqn5 = diff(Ppv(t),t) == (((Ppp(t) - Ppv(t))/Rpp)...
         - ((Ppv(t) - Pla(t))/Rpv))/Cpv;
eqn6 = diff(Fsa(t)) == (Psa(t) - Psp(t) - Rsa*Fsa(t));
eqn7 = diff(Psp(t),t) == (1/(Csp+Cep)) * (Fsa(t) - ((Psp(t) ...
    - Psv(t))/Rsp) - ((Psp(t) - Pev(t))/Rep));
eqn8 = diff(Pev(t),t) == (1/Cev) * (((Psp(t) - Pev(t))/Rep)...
    - ((Pev(t) - Pra(t))/Rev) - diff(Vuev(t),t));
eqn9 = Psv(t) == (1/Csv) * (Vt(t) - Csa*Psa(t) - (Csp+Cep)*Psp(t)...
    - Cev*Pev(t) - Cra*Pra(t) - Vrv(t) - Cpa*Ppa(t) - Cpp*Ppp(t)...
    - Cpv*Ppv(t) - Cla*Pla(t) - Vlv(t) - Vu(t));
eqn10 = Vu(t) == Vusa + Vusp + Vuep + Vusv + Vuev + Vura + Vupa + Vupp...
    + Vupv + Vula;
eqn11 = diff(Vt(t),t) == I(t);
eqn12 = Vt(0) == Vt0;
eqn13 = diff(Pla(t),t) == (1/Cla) * (((Ppv(t) - Pla(t))/Rpv) - Fil);

if (Pla(t) <= Plv(t))
    eqn14 = Fil(t) == 0;
else
    eqn14 = Fil(t) == ((Pla(t) - Plv(t))/Rla);
end

eqn15 = diff(Vlv(t),t) == Fil(t) - Fol(t);

if (Pmaxlv(t) <= Psa(t))
    eqn16 = Fol(t) == 0;
else
    eqn16 = Fol(t) == ((Pmaxlv(t) - Psa(t))/Rlv(t));
end

eqn17 = Rlv(t) == kRlv*Pmaxlv(t);
eqn18 = Plv(t) == Pmaxlv(t) - Rlv(t)*Fol(t);
eqn19 = Pmaxlv(t) == phi(t) * Emaxlv * (Vlv(t) - Vulv)...
    + (1 - phi(t))*P0lv*(exp(kElv*Vlv(t)) - 1);

if (0 <= u(t) && u(t) <= (Tsys(t) / T(t)))
    eqn20 = phi(t) == (sin(((pi * T(t))/Tsys(t)) * u(t)))^2;
else
    eqn20 = phi(t) == 0;
end

eqn21 = u(t) == mod(abs(fix((int( (1/T(t)), [t0, t] ) + u(t0)))), 1);
eqn22 = diff(epsilon,t) == (1/T(t));
eqn23 = u(t) == mod(abs(epsilon), 1);
eqn24 = Tsys(t) == Tsys0 - ksys*(1/T(t));
eqn25 = diff(Pra(t),t) == (1/Cra) * (((Psv(t) - Pra(t))/Rsv)...
        + ((Pev(t) - Pra(t))/Rev) - Fir(t));

if (Pra(t) <= Prv(t))
    eqn26 = Fir(t) == 0;
else
    eqn26 = Fir(t) == ((Pra(t) - Prv(t))/Rra);
end

eqn27 = diff(Vrv(t),t) == Fir(t) - For(t);

if (Pmaxlv(t) <= Psa(t))
    eqn28 = For(t) == 0;
else
    eqn28 = For(t) == ((Pmaxrv(t) - Ppa(t))/Rrv(t));
end

eqn29 = Rrv(t) == kRrv * Pmaxrv(t);
eqn30 = Prv(t) == Pmaxrv(t) - Rrv(t)*For(t);
eqn31 = Pmaxrv(t) == phi(t) * Emaxrv * (Vrv(t) - Vurv)...
    + (1 - phi(t))*P0rv*(exp(kErv*Vrv(t)) - 1);
eqn32 = taup*diff(Ptil(t),t) == Pcs(t) + tauz*diff(Pcs(t),t) - Ptil(t);
eqn33 = fcs(t) == (fmin + fmax*exp((Ptil(t) - Pn)/ka))...
    / (1 + exp((Ptil(t) - Pn)/ka));
eqn34 = ka == ((fmax - fmin) / (4*Gb));
eqn35 = fes(t) == fes_8 + (fes0 - fes_8)*exp(-kes*fcs(t));
eqn36 = fev(t) == (fev0 + fev_8*exp((fcs(t) - fcs0)/kev))...
    / (1 + exp((fcs(t) - fcs0)/kev));

if (fes(t) >= fesmin)
    eqn37 = sigma_Emaxrv(t) == GEmaxrv * log(fes(t - DEmaxrv) - fesmin...
    +1);
else
    eqn37 = sigma_Emaxrv(t) == 0;
end

if (fes(t) >= fesmin)
    eqn38 = sigma_Emaxlv(t) == GEmaxlv * log(fes(t - DEmaxlv) - fesmin...
    +1);
else
    eqn38 = sigma_Emaxlv(t) == 0;
end

if (fes(t) >= fesmin)
    eqn39 = sigma_Rsp(t) == GRsp * log(fes(t - DRsp) - fesmin...
    +1);
else
    eqn39 = sigma_Rsp(t) == 0;
end

if (fes(t) >= fesmin)
    eqn40 = sigma_Rep(t) == GRep * log(fes(t - DRep) - fesmin...
    +1);
else
    eqn40 = sigma_Rep(t) == 0;
end

if (fes(t) >= fesmin)
    eqn41 = sigma_Vusv(t) == GVusv * log(fes(t - DVusv) - fesmin...
    +1);
else
    eqn41 = sigma_Vusv(t) == 0;
end

if (fes(t) >= fesmin)
    eqn42 = sigma_Vuev(t) == GVuev * log(fes(t - DVuev) - fesmin...
    +1);
else
    eqn42 = sigma_Vuev(t) == 0;
end

eqn43 = diff(delta_Emaxrv(t),t) == (1/tauEmaxrv) * (-delta_Emaxrv(t)...
    + sigma_Emaxrv(t));
eqn44 = diff(delta_Emaxlv(t),t) == (1/tauEmaxlv) * (-delta_Emaxlv(t)...
    + sigma_Emaxlv(t));
eqn45 = diff(delta_Rsp(t),t) == (1/tauRsp) * (-delta_Rsp(t)...
    + sigma_Rsp(t));
eqn46 = diff(delta_Rep(t),t) == (1/tauRep) * (-delta_Rep(t)...
    + sigma_Rep(t));
eqn47 = diff(delta_Vusv(t),t) == (1/tauVusv) * (-delta_Vusv(t)...
    + sigma_Vusv(t));
eqn49 = diff(delta_Vuev(t),t) == (1/tauVuev) * (-delta_Vuev(t)...
    + sigma_Vuev(t));

eqn50 = Emaxrv == delta_Emaxrv(t) + Emaxrv0;
eqn51 = Emaxlv == delta_Emaxlv(t) + Emaxlv0;
eqn52 = Rsp == delta_Rsp(t) + Rsp0;
eqn53 = Rep == delta_Rep(t) + Rep0;
eqn54 = Vusv == delta_Vusv(t) + Vusv0;
eqn55 = Vuev == delta_Vuev(t) + Vuev0;

if (fes(t) >= fesmin)
    eqn56 = sigma_Ts(t) == GTs * log(fes(t - DTs) - fesmin + 1);
else
    eqn56 = sigma_Ts(t) == 0;
end

eqn57 = diff(delta_Ts(t),t) == (1/tauTs) * (-delta_Ts(t) + sigma_Ts(t));
eqn58 = sigma_Tv(t) == GTv * fev(t - DTv);
eqn59 = diff(delta_Tv(t),t) == (1/tauTv) * (-delta_Tv(t) + sigma_Tv(t));
eqn60 = T(t) == delta_Ts(t) + delta_Tv(t) + T0;

end