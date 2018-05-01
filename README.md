# PairPurityFitter



```
uls_pid_m27 = uls_pid_m27
template_uls_sumbg_m27 = template_uls_sumbg_m27
ls_pid_m27 = ls_pid_m27
template_ls_sumbg_m27 = template_ls_sumbg_m27

template_uls_sumbg_m27->Scale( uls_pid_m27->Integral() )
template_uls_sumbg_m27->Draw("same")

template_ls_sumbg_m27->Scale( ls_pid_m27->Integral() )
template_ls_sumbg_m27->Draw("same")
template_ls_sumbg_m27->SetLineColor(kRed)
template_uls_sumbg_m27->Divide( template_ls_sumbg_m27 )

uls_pid_m27->Divide( ls_pid_m27 )
uls_pid_m27->Draw("h")

template_uls_sumbg_m27->Draw("same hist")


```

