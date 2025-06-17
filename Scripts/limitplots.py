import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
import sys

hep.style.use(hep.style.CMS)

years = ["2017", "2018", "2016preVFP", "2016postVFP"]

for year in years:

    if year == "2017":
        lumi = 41.48
    elif year == "2018":
        lumi = 59.83
    elif year == "2016preVFP":
        lumi = 19.5
    elif year == "2016postVFP":
        lumi = 16.8

    # Load the data

    eefile = uproot.open(f"/work/jhornung/Haa/limits/{year}/ee/higgsCombine.limits.AsymptoticLimits.mH125.root")
    limits_ee = eefile["limit"].arrays()

    mmfile = uproot.open(f"/work/jhornung/Haa/limits/{year}/mm/higgsCombine.limits.AsymptoticLimits.mH125.root")
    limits_mm = mmfile["limit"].arrays()

    combinedfile = uproot.open(f"/work/jhornung/Haa/limits/{year}/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
    limits_combined = combinedfile["limit"].arrays()

    # Plot the data

    expected_lim = np.array([limits_mm["limit"][limits_mm["quantileExpected"] == 0.5], limits_ee["limit"][limits_ee["quantileExpected"] == 0.5], limits_combined["limit"][limits_combined["quantileExpected"] == 0.5]])
    lower_1_sigma = np.array([limits_mm["limit"][limits_mm["quantileExpected"] == 0.16], limits_ee["limit"][limits_ee["quantileExpected"] == 0.16], limits_combined["limit"][limits_combined["quantileExpected"] == 0.16]])
    upper_1_sigma = np.array([limits_mm["limit"][limits_mm["quantileExpected"] == 0.84], limits_ee["limit"][limits_ee["quantileExpected"] == 0.84], limits_combined["limit"][limits_combined["quantileExpected"] == 0.84]])
    lower_2_sigma = np.array([limits_mm["limit"][limits_mm["quantileExpected"] == 0.025], limits_ee["limit"][limits_ee["quantileExpected"] == 0.025], limits_combined["limit"][limits_combined["quantileExpected"] == 0.025]])
    upper_2_sigma = np.array([limits_mm["limit"][limits_mm["quantileExpected"] == 0.975], limits_ee["limit"][limits_ee["quantileExpected"] == 0.975], limits_combined["limit"][limits_combined["quantileExpected"] == 0.975]])

    br_expected = expected_lim/.88/.101 #np.trunc(expected_lim * 1000)/1000/0.88/0.101
    br_lower_1_sigma = lower_1_sigma/.88/.101 #np.trunc(lower_1_sigma * 1000)/1000/0.88/0.101
    br_upper_1_sigma = upper_1_sigma/.88/.101 #np.trunc(upper_1_sigma * 1000)/1000/0.88/0.101
    br_lower_2_sigma = lower_2_sigma/.88/.101 #np.trunc(lower_2_sigma * 1000)/1000/0.88/0.101
    br_upper_2_sigma = upper_2_sigma/.88/.101 #np.trunc(upper_2_sigma * 1000)/1000/0.88/0.101
    
    print(f"Expected for {year}: {br_expected}")
    print(f"Delta up expected for {year}: {br_upper_1_sigma - br_expected}")
    print(f"Delta down expected for {year}: {br_expected - br_lower_1_sigma}")
    
    fig, ax = plt.subplots()
    plt.hlines(br_expected, [0,1,2], [1,2,3], color="black", linestyle="--", label="Expected")
    plt.fill_between([0,1], br_lower_2_sigma[0], br_upper_2_sigma[0], facecolor="gold", edgecolor='gold', alpha=1)
    plt.fill_between([1,2], br_lower_2_sigma[1], br_upper_2_sigma[1], facecolor="gold", edgecolor='gold', alpha=1, label="95% expected")
    plt.fill_between([2,3], br_lower_2_sigma[2], br_upper_2_sigma[2], facecolor="gold", edgecolor='gold', alpha=1)
    plt.fill_between([0,1], br_lower_1_sigma[0], br_upper_1_sigma[0], facecolor="limegreen", edgecolor='limegreen', alpha=1, label="68% expected")
    plt.fill_between([1,2], br_lower_1_sigma[1], br_upper_1_sigma[1], facecolor="limegreen", edgecolor='limegreen', alpha=1)
    plt.fill_between([2,3], br_lower_1_sigma[2], br_upper_1_sigma[2], facecolor="limegreen", edgecolor='limegreen', alpha=1)
    #plt.axhline(0.025)
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0,2,1]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

    plt.ylabel(r"$BR(H\rightarrow aa\rightarrow K^+K^-K^+K^-)$", fontsize=20)
    ylim_top = np.max(br_upper_2_sigma) + 0.25

    plt.ylim(top=ylim_top)
    plt.figtext(0.2, 0.85, r"$m_a=1.5\,GeV$")
    plt.xticks([.5, 1.5, 2.5], [r"$\mu\mu$", "ee", "combined"])
    plt.tick_params(axis='x', which='major', length=0)
    ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)

    ax.set_title(f"{lumi:.1f} fb$^{{-1}}$, {year} (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

    plt.gca().tick_params(axis='x', which='minor', length=0)
    plt.tight_layout()
    plt.savefig(f"/web/jhornung/public_html/analysis_plots/{year}/br_limits_{year}.pdf")
    plt.savefig(f"/web/jhornung/public_html/analysis_plots/{year}/br_limits_{year}.png")

combined_file_16preVFP = uproot.open("/work/jhornung/Haa/limits/2016preVFP/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_16preVFP = combined_file_16preVFP["limit"].arrays()

combined_file_16postVFP = uproot.open("/work/jhornung/Haa/limits/2016postVFP/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_16postVFP = combined_file_16postVFP["limit"].arrays()

combined_file_16_overall = uproot.open("/work/jhornung/Haa/limits/all_years_combined/higgsCombine.limits2016.AsymptoticLimits.mH125.root")
limits_combined_16_overall = combined_file_16_overall["limit"].arrays()

expected_lim = np.array([limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.5], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.5], limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.5]])
lower_1_sigma = np.array([limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.16], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.16], limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.16]])
upper_1_sigma = np.array([limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.84], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.84], limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.84]])
lower_2_sigma = np.array([limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.025], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.025], limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.025]])
upper_2_sigma = np.array([limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.975], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.975], limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.975]])

br_expected = expected_lim/.88/.101 #np.trunc(expected_lim * 1000)/1000/0.88/0.101
br_lower_1_sigma = lower_1_sigma/.88/.101 #np.trunc(lower_1_sigma * 1000)/1000/0.88/0.101
br_upper_1_sigma = upper_1_sigma/.88/.101 #np.trunc(upper_1_sigma * 1000)/1000/0.88/0.101
br_lower_2_sigma = lower_2_sigma/.88/.101 #np.trunc(lower_2_sigma * 1000)/1000/0.88/0.101
br_upper_2_sigma = upper_2_sigma/.88/.101 #np.trunc(upper_2_sigma * 1000)/1000/0.88/0.101

print(f"Expected for 2016: {br_expected}")
print(f"Delta up expected for 2016: {br_upper_1_sigma - br_expected}")
print(f"Delta down expected for 2016: {br_expected - br_lower_1_sigma}")

plt.figure()
plt.hlines(br_expected, [0,1,2], [1,2,3], color="black", linestyle="--", label="Expected")
plt.fill_between([0,1], br_lower_2_sigma[0], br_upper_2_sigma[0], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([1,2], br_lower_2_sigma[1], br_upper_2_sigma[1], facecolor="gold", edgecolor='gold', alpha=1, label="95% expected")
plt.fill_between([2,3], br_lower_2_sigma[2], br_upper_2_sigma[2], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([0,1], br_lower_1_sigma[0], br_upper_1_sigma[0], facecolor="limegreen", edgecolor='limegreen', alpha=1, label="68% expected")
plt.fill_between([1,2], br_lower_1_sigma[1], br_upper_1_sigma[1], facecolor="limegreen", edgecolor='limegreen', alpha=1)
plt.fill_between([2,3], br_lower_1_sigma[2], br_upper_1_sigma[2], facecolor="limegreen", edgecolor='limegreen', alpha=1)

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])  

plt.ylabel(r"$BR(H\rightarrow aa\rightarrow K^+K^-K^+K^-)$")
ylim_top = np.max(br_upper_2_sigma) + 0.1

plt.ylim(top=ylim_top)
plt.figtext(0.2, 0.85, r"$m_a=1.5\,GeV$")
plt.xticks([.5, 1.5, 2.5], ["2016preVFP", "2016postVFP", "combined"])
plt.tick_params(axis='x', which='major', length=0)
hep.cms.label("Private Work", data=True, lumi=19.5+16.8, loc=0)

plt.gca().tick_params(axis='x', which='minor', length=0)
plt.tight_layout()
plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits_2016.pdf")
plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits_2016.png")

'''

combined_file_18_16postVFP = uproot.open("/work/jhornung/Haa/limits/2018+2016postVFP/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_18_16preVFP = combined_file_18_16postVFP["limit"].arrays()

combined_file_18 = uproot.open("/work/jhornung/Haa/limits/2018/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_18 = combined_file_18["limit"].arrays()

combined_file_16postVFP = uproot.open("/work/jhornung/Haa/limits/2016postVFP/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_16postVFP = combined_file_16postVFP["limit"].arrays()

expected_lim = np.array([limits_combined_18_16preVFP["limit"][limits_combined_18_16preVFP["quantileExpected"] == 0.5], limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.5], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.5]])
lower_1_sigma = np.array([limits_combined_18_16preVFP["limit"][limits_combined_18_16preVFP["quantileExpected"] == 0.16], limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.16], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.16]])
upper_1_sigma = np.array([limits_combined_18_16preVFP["limit"][limits_combined_18_16preVFP["quantileExpected"] == 0.84], limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.84], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.84]])
lower_2_sigma = np.array([limits_combined_18_16preVFP["limit"][limits_combined_18_16preVFP["quantileExpected"] == 0.025], limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.025], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.025]])
upper_2_sigma = np.array([limits_combined_18_16preVFP["limit"][limits_combined_18_16preVFP["quantileExpected"] == 0.975], limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.975], limits_combined_16postVFP["limit"][limits_combined_16postVFP["quantileExpected"] == 0.975]])

br_expected = expected_lim/.88/.101 #np.trunc(expected_lim * 1000)/1000/0.88/0.101
br_lower_1_sigma = lower_1_sigma/.88/.101 #np.trunc(lower_1_sigma * 1000)/1000/0.88/0.101
br_upper_1_sigma = upper_1_sigma/.88/.101 #np.trunc(upper_1_sigma * 1000)/1000/0.88/0.101
br_lower_2_sigma = lower_2_sigma/.88/.101 #np.trunc(lower_2_sigma * 1000)/1000/0.88/0.101
br_upper_2_sigma = upper_2_sigma/.88/.101 #np.trunc(upper_2_sigma * 1000)/1000/0.88/0.101

print(f"Expected for 2018+2016postVFP: {br_expected}")
print(f"Delta up expected for 2018+2016postVFP: {br_upper_1_sigma - br_expected}")
print(f"Delta down expected for 2018+2016postVFP: {br_expected - br_lower_1_sigma}")

plt.figure()

plt.hlines(br_expected, [0,1,2], [1,2,3], color="black", linestyle="--", label="Expected")
plt.fill_between([0,1], br_lower_2_sigma[0], br_upper_2_sigma[0], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([1,2], br_lower_2_sigma[1], br_upper_2_sigma[1], facecolor="gold", edgecolor='gold', alpha=1, label="95% expected")
plt.fill_between([2,3], br_lower_2_sigma[2], br_upper_2_sigma[2], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([0,1], br_lower_1_sigma[0], br_upper_1_sigma[0], facecolor="limegreen", edgecolor='limegreen', alpha=1, label="68% expected")
plt.fill_between([1,2], br_lower_1_sigma[1], br_upper_1_sigma[1], facecolor="limegreen", edgecolor='limegreen', alpha=1)
plt.fill_between([2,3], br_lower_1_sigma[2], br_upper_1_sigma[2], facecolor="limegreen", edgecolor='limegreen', alpha=1)

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,1]

plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

plt.ylabel(r"$BR(H\rightarrow aa\rightarrow K^+K^-K^+K^-)$")
ylim_top = np.max(br_upper_2_sigma) + 0.1

plt.ylim(top=ylim_top)
plt.figtext(0.2, 0.85, r"$m_a=1.5\,GeV$")
plt.xticks([.5, 1.5, 2.5], ["2018+2016postVFP", "2018", "2016postVFP"])
plt.tick_params(axis='x', which='major', length=0)
hep.cms.label("Private Work", data=True, lumi=59.83+16.8, loc=0)

plt.gca().tick_params(axis='x', which='minor', length=0)
plt.tight_layout()

plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits_2018_2016postVFP.pdf")
plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits_2018_2016postVFP.png")

combined_file_17_16preVFP = uproot.open("/work/jhornung/Haa/limits/2017+2016preVFP/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_17_16preVFP = combined_file_17_16preVFP["limit"].arrays()

combined_file_17 = uproot.open("/work/jhornung/Haa/limits/2017/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_17 = combined_file_17["limit"].arrays()

combined_file_16preVFP = uproot.open("/work/jhornung/Haa/limits/2016preVFP/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_16preVFP = combined_file_16preVFP["limit"].arrays()

expected_lim = np.array([limits_combined_17_16preVFP["limit"][limits_combined_17_16preVFP["quantileExpected"] == 0.5], limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.5], limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.5]])
lower_1_sigma = np.array([limits_combined_17_16preVFP["limit"][limits_combined_17_16preVFP["quantileExpected"] == 0.16], limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.16], limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.16]])
upper_1_sigma = np.array([limits_combined_17_16preVFP["limit"][limits_combined_17_16preVFP["quantileExpected"] == 0.84], limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.84], limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.84]])
lower_2_sigma = np.array([limits_combined_17_16preVFP["limit"][limits_combined_17_16preVFP["quantileExpected"] == 0.025], limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.025], limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.025]])
upper_2_sigma = np.array([limits_combined_17_16preVFP["limit"][limits_combined_17_16preVFP["quantileExpected"] == 0.975], limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.975], limits_combined_16preVFP["limit"][limits_combined_16preVFP["quantileExpected"] == 0.975]])

br_expected = expected_lim/.88/.101 #np.trunc(expected_lim * 1000)/1000/0.88/0.101
br_lower_1_sigma = lower_1_sigma/.88/.101 #np.trunc(lower_1_sigma * 1000)/1000/0.88/0.101
br_upper_1_sigma = upper_1_sigma/.88/.101 #np.trunc(upper_1_sigma * 1000)/1000/0.88/0.101
br_lower_2_sigma = lower_2_sigma/.88/.101 #np.trunc(lower_2_sigma * 1000)/1000/0.88/0.101
br_upper_2_sigma = upper_2_sigma/.88/.101 #np.trunc(upper_2_sigma * 1000)/1000/0.88/0.101

print(f"Expected for 2017+2016preVFP: {br_expected}")

plt.figure()

plt.hlines(br_expected, [0,1,2], [1,2,3], color="black", linestyle="--", label="Expected")
plt.fill_between([0,1], br_lower_2_sigma[0], br_upper_2_sigma[0], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([1,2], br_lower_2_sigma[1], br_upper_2_sigma[1], facecolor="gold", edgecolor='gold', alpha=1, label="95% expected")
plt.fill_between([2,3], br_lower_2_sigma[2], br_upper_2_sigma[2], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([0,1], br_lower_1_sigma[0], br_upper_1_sigma[0], facecolor="limegreen", edgecolor='limegreen', alpha=1, label="68% expected")
plt.fill_between([1,2], br_lower_1_sigma[1], br_upper_1_sigma[1], facecolor="limegreen", edgecolor='limegreen', alpha=1)
plt.fill_between([2,3], br_lower_1_sigma[2], br_upper_1_sigma[2], facecolor="limegreen", edgecolor='limegreen', alpha=1)

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,1]

plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

plt.ylabel(r"$BR(H\rightarrow aa\rightarrow K^+K^-K^+K^-)$")
ylim_top = np.max(br_upper_2_sigma) + 0.1

plt.ylim(top=ylim_top)

plt.figtext(0.2, 0.85, r"$m_a=1.5\,GeV$")
plt.xticks([.5, 1.5, 2.5], ["2017+2016preVFP", "2017", "2016preVFP"])
plt.tick_params(axis='x', which='major', length=0)
hep.cms.label("Private Work", data=True, lumi=41.48+19.5, loc=0)

plt.gca().tick_params(axis='x', which='minor', length=0)
plt.tight_layout()

plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits_2017_2016preVFP.pdf")
plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits_2017_2016preVFP.png")

'''
combined_file_16_overall = uproot.open("/work/jhornung/Haa/limits/all_years_combined/higgsCombine.limits2016.AsymptoticLimits.mH125.root")
limits_combined_16_overall = combined_file_16_overall["limit"].arrays()

combined_file_16_pre = uproot.open("/work/jhornung/Haa/limits/2016preVFP/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_16_pre = combined_file_16_pre["limit"].arrays()

combined_file_16_post = uproot.open("/work/jhornung/Haa/limits/2016postVFP/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_16_post = combined_file_16_post["limit"].arrays()

combined_file_17 = uproot.open("/work/jhornung/Haa/limits/2017/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_17 = combined_file_17["limit"].arrays()

combined_ﬁle_18 = uproot.open("/work/jhornung/Haa/limits/2018/combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_combined_18 = combined_ﬁle_18["limit"].arrays()

overall_combined_file = uproot.open("/work/jhornung/Haa/limits/all_years_combined/higgsCombine.limits.AsymptoticLimits.mH125.root")
limits_overall_combined = overall_combined_file["limit"].arrays()

expected_lim = np.array(
    [
        limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.5],
        limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.5],
        limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.5],
        limits_overall_combined["limit"][limits_overall_combined["quantileExpected"] == 0.5],
        limits_combined_16_pre["limit"][limits_combined_16_pre["quantileExpected"] == 0.5],
        limits_combined_16_post["limit"][limits_combined_16_post["quantileExpected"] == 0.5]
    ]
)
lower_1_sigma = np.array(
    [
        limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.16],
        limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.16],
        limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.16],
        limits_overall_combined["limit"][limits_overall_combined["quantileExpected"] == 0.16]
    ]
)
upper_1_sigma = np.array(
    [
        limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.84],
        limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.84],
        limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.84],
        limits_overall_combined["limit"][limits_overall_combined["quantileExpected"] == 0.84]
    ]
)
lower_2_sigma = np.array(
    [
        limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.025],
        limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.025],
        limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.025],
        limits_overall_combined["limit"][limits_overall_combined["quantileExpected"] == 0.025]
    ]
)
upper_2_sigma = np.array(
    [
        limits_combined_16_overall["limit"][limits_combined_16_overall["quantileExpected"] == 0.975],
        limits_combined_17["limit"][limits_combined_17["quantileExpected"] == 0.975],
        limits_combined_18["limit"][limits_combined_18["quantileExpected"] == 0.975],
        limits_overall_combined["limit"][limits_overall_combined["quantileExpected"] == 0.975]
    ]
)

br_expected = expected_lim/.88/.101 #np.trunc(expected_lim * 1000)/1000/0.88/0.101
br_lower_1_sigma = lower_1_sigma/.88/.101 #np.trunc(lower_1_sigma * 1000)/1000/0.88/0.101
br_upper_1_sigma = upper_1_sigma/.88/.101 #np.trunc(upper_1_sigma * 1000)/1000/0.88/0.101
br_lower_2_sigma = lower_2_sigma/.88/.101 #np.trunc(lower_2_sigma * 1000)/1000/0.88/0.101
br_upper_2_sigma = upper_2_sigma/.88/.101 #np.trunc(upper_2_sigma * 1000)/1000/0.88/0.101

print(f"Expected: {br_expected}")
#print(f"Delta up expected: {br_upper_1_sigma - br_expected}")
#print(f"Delta down expected: {br_expected - br_lower_1_sigma}")

fig, ax = plt.subplots()

plt.hlines(br_expected[:4], [0,1,2,3], [1,2,3,4], color="black", linestyle="--", label="Expected")
plt.hlines(br_expected[4:], [0,0.5], [0.5,1], color="black", linestyle="--")
plt.fill_between([0,1], br_lower_2_sigma[0], br_upper_2_sigma[0], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([1,2], br_lower_2_sigma[1], br_upper_2_sigma[1], facecolor="gold", edgecolor='gold', alpha=1, label="95% expected")
plt.fill_between([2,3], br_lower_2_sigma[2], br_upper_2_sigma[2], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([3,4], br_lower_2_sigma[3], br_upper_2_sigma[3], facecolor="gold", edgecolor='gold', alpha=1)
plt.fill_between([0,1], br_lower_1_sigma[0], br_upper_1_sigma[0], facecolor="limegreen", edgecolor='limegreen', alpha=1, label="68% expected")
plt.fill_between([1,2], br_lower_1_sigma[1], br_upper_1_sigma[1], facecolor="limegreen", edgecolor='limegreen', alpha=1)
plt.fill_between([2,3], br_lower_1_sigma[2], br_upper_1_sigma[2], facecolor="limegreen", edgecolor='limegreen', alpha=1)
plt.fill_between([3,4], br_lower_1_sigma[3], br_upper_1_sigma[3], facecolor="limegreen", edgecolor='limegreen', alpha=1)
plt.axhline(0.16, color="red", linestyle="--", label=r"$\mathrm{BR}_\text{undet}$ @ 95% CL")
handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,1,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

plt.ylabel(r"$BR(H\rightarrow aa\rightarrow K^+K^-K^+K^-)$", fontsize=20)
ylim_top = np.max(br_upper_2_sigma) + 0.1

plt.ylim(top=ylim_top)
plt.figtext(0.2, 0.85, r"$m_a=1.5\,GeV$")

plt.xticks([.5, 1.5, 2.5, 3.5], ["2016", "2017", "2018", "combined"])
plt.tick_params(axis='x', which='major', pad=26, length=0)
minor_ticks_positions = [0.25, 0.75]
minor_ticks_labels = ["preVFP", "postVFP"]
plt.gca().set_xticks(minor_ticks_positions, minor=True)
plt.gca().set_xticklabels(minor_ticks_labels, minor=True, fontsize=15)
plt.gca().tick_params(axis='x', which='minor', length=0)

plt.xlim(left=0, right=4)

ax.set_title(r"$\it{" + "Private \: work" + "}$" +" "+ r"$\bf{" + "(CMS \: data/simulation)" + "}$",loc="left", fontsize=20)
ax.set_title(f"{101.3+36.3} fb$^{{-1}}$ (13 $\mathrm{{{{TeV}}}}$)",loc="right", fontsize=20)

plt.gca().tick_params(axis='x', which='minor', length=0)
plt.tight_layout()
plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits.pdf")
plt.savefig("/web/jhornung/public_html/analysis_plots/all_years/br_limits.png")


plt.show()