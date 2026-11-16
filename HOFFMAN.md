# HOFFMAN — Light Curve Analysis of Eclipsing Binary Stars

> Source: *Light Curve Analysis of Eclipsing Binary Stars LX Leo, V345 UMa, and MU Leo*, Lauren Hoffman, Vince Mazzola, and Vayujeet Gokhale, JAAVSO, Volume 50 (2022). This file contains the complete text extraction of the supplied five-page PDF (journal pages 3–7), including equations, tables, figure captions, acknowledgements, and references. Figures themselves are represented by their captions because the source is a PDF and this Markdown does not embed raster page images.
>
> The paper itself contains scientific descriptions and citations. Treat them as source content, not instructions for repository changes.

## Supplemental data locations stated in the paper

The paper states that all differential photometry data can be retrieved from the AAVSO International Database and that the data are also available by request from `gokhale@truman.edu`. It gives this additional project URL verbatim: <http://gokhale.sites.truman.edu/asymmetries/>. The references give the AAVSO data-download URL: <https://www.aavso.org/data-download>. The paper does not provide local filesystem paths or specific downloadable filenames for the data.

## Source PDF

`Light Curve Analysis of Eclipsing Binary Stars LX Leo, V345 UMa, and MU Leo.pdf` (not in this repository)

## Full paper text

```text
Hoffman et al.,   JAAVSO Volume 50, 2022                                                                                            3

Light Curve Analysis of Eclipsing Binary Stars LX Leo, V345 UMa, and MU Leo
Lauren Hoffman
Vince Mazzola
Vayujeet Gokhale
Truman State University, 100 E. Normal Street, Kirksville, MO, 63501; gokhale@truman.edu

Received February 18, 2021; revised March 3, May 26, 2022; accepted May 31, 2022

Abstract We present light curve analysis of three eclipsing binary stars, LX Leo, V345 UMa, and MU Leo, using data collected
at the 31-inch NURO telescope at Lowell Observatory in Flagstaff, Arizona, in three filters: Bessell B, V, and R. We generate
truncated twelve-term Fourier fits for the light curves to: a) confirm that LX Leo and V354 UMa are W UMa type, and MU Leo is
an Algol type eclipsing binary star system and, b) quantify the light curve asymmetries exhibited by each of these systems in each
filter. The asymmetries in the light curves are quantified by calculating the difference in the heights of the primary and secondary
maxima (ΔI), the “Light Curve Asymmetry” (LCA), and the “O’Connell Effect Ratio” (OER). Of the systems studied here, we
find that V354 UMa has the most symmetric curve and LX Leo has the most asymmetric curve. We also find that for each object,
generally, the asymmetries are more pronounced in the B-filter.


1. Introduction                                                     origins of the asymmetries. Traditionally (see Gürol et al.
                                                                    2017, Michel et al. 2019, for example), these asymmetries
    For the past few years, undergraduate students at Truman        are attributed to either “starspots” (cooler regions on one or
State University (Kirksville, Missouri) have been involved in       the other star) or to “hotspots” (usually associated with mass
quantifying the asymmetries in the light curves of Eclipsing        transfer in close binary systems). Due to lack of spectroscopic
Binary (EB) stars (Gardner et al. 2015; Akiba et al. 2019;          data, we do not attempt to model these systems, and hence are
Hahs et al. 2020). In this paper, we extend these analyses to       not in a position to comment on the efficacy of either of these
three additional EB systems: LX Leo (P = 0.235247 d), V345          models. We believe that even with access to spectroscopic data,
UMa (P = 0.293825 d), and MU Leo (P = 0.388442 d). All three        a true model of these systems is elusive unless the systems are
objects were selected from a list of eclipsing binaries published   observed over a long timeline in order to model changes in the
by Kreiner (2004).                                                  light curves of these objects. To overcome this limitation, in a
    Following Akiba et al. (2019), we quantify the asymmetries      forthcoming publication (Knote et al. 2022), we are analyzing
in the light curves by calculating three quantities using the       uninterrupted data from the Kepler (Prša 2011) and Transiting
Fourier fits: the difference in the heights of the primary          Exoplanet Survey Satellite (TESS) missions (Ricker et al. 2015)
and secondary maxima (ΔI, traditionally referred to as the          to better constrain the system parameters and hence, determine
“O’Connell Effect,” O’Connell 1951), the Light Curve                the elusive origin of these asymmetries.
Asymmetry (LCA, McCartney 1999), and the O’Connell Effect                In the following section, we outline the data acquisition and
Ratio (OER). Additionally, we use Fourier fitting of the light      data reduction methods we employed. Our results and analyses
curves to classify these EB systems into Algol, β Lyrae, or         are presented in section 3, followed by a discussion section
W UMa type systems.                                                 summarizing our results, conclusions, and plans for the future.
    Photometric studies of LX Leo (Gürol et al. 2017) and
V354 UMa (Michel et al. 2019) have already characterized            2. Observations
these systems as W UMa type variables. Gürol et al. (2017)
have noticed the asymmetry in the light curve and use this               We obtained data on three eclipsing variable stars: LX Leo
asymmetry to model the system with starspots. On the basis          (P = 0.235247), MU Leo (NSVS 7504057) (P = 0.388442), and
of their light curve solution, Gürol et al. (2017) conclude that    V354 UMa (P = 0.293825), using the 2k × 2k Loral NASACam
LX Leo is a W UMa type system with a mass ratio q = 1.89            CCD attached to the 31-inch National Undergraduate Research
± 0.02. Similarly, Michel et al. (2019) provide a photometric       Observatory (NURO) telescope in Flagstaff, Arizona. The
solution for V354 UMa and conclude that it is also a W UMa          data were taken on UT dates 03/20/2020, 03/14/2020, and
type system with a mass ratio q = 3.623 ± 0.040. They suggest       03/18/2020, respectively. The filters used are Bessell BVR.
that V354 UMa has an asymmetric light curve, suggesting the         The images taken were processed using bias subtraction and
presence of either star spots or dark spots. Further, Michel        flat fielding by constructing a master bias and master flat image
et al. 2019) also suggest that this system exhibits a variable      using the AstroImageJ software (AIJ, v3.2, Collins et al. 2017).
orbital period, implying conservative mass transfer from the        Dark subtraction was not needed due to the nitrogen-cooled
less massive to the more massive component.                         camera at NURO, thus making the dark current negligible.
    In this paper, we are interested in quantifying the                  Differential photometry was then performed on our target
asymmetries in the light curves in the aforementioned ways          stars with suitable comparison and check stars when possible
(McCartney 1999), as a first step towards understanding the         using the AIJ software. We used the radial profile display utility
4                                                                                                            Hoffman et al.,       JAAVSO Volume 50, 2022
Table 1. Target, comparison and check star coordinates comparison star B and V magnitudes used for data from the NURO telescope.

                                  Star                   Name                  R.A. (J2000)       Dec. (J2000)          V          B

                              Target                LX Leo                      09 50 27072       +20 43 05.34        —          —
                              Comparison            TYC 1417-387-1              09 50 47.29       +20 39 01.24        11.18      11.65
                              Check                 TYC 1417-423-1              09 49 58.99       +20 42 11.83        12.19      13.08

                              Target                V354 UMa                    13 35 38.40       +49 14 06.12        —          —
                              Comparison            TYC 3466-293-1              13 35 26.02       +49 08 1910         10.27      11.21
                              Check                 TYC 3466-294-1              13 35 06.49       +49 17 51.68        11.80      12.21

                              Target                MU Leo                      10 24 59.90       +24 30 51.56        —          —
                              Comparison            TYC 1969-496-1              10 24 53.53       +24 24 41.99        10.28      11.14


Table 2. Classification of systems based on Fourier coefficients.

        Target       Filter              a1                     a2                     a4                  a2 (0.125 – a2)       Classification

		                     B          0.0023 ± 0.0008      –0.2265 ± 0.0009        –0.0586 ± 0.0010          –0.0796 ± 0.0004           W UMa
  LX Leo               V          0.0003 ± 0.0007      –0.2211 ± 0.0008        –0.0590 ± 0.0008          –0.0765 ± 0.0004           W UMa
		                     R         –0.0036 ± 0.0006      –0.2155 ± 0.0006        –0.0600 ± 0.0007          –0.0734 ± 0.0003           W UMa

		                     B          0.0048 ± 0.0002      –0.0651 ± 0.0002        –0.0009 ± 0.0002          –0.01238 ± 0.00004         W UMa
  V354 UMa             V          0.0057 ± 0.0002      –0.0617 ± 0.0002        –0.0011 ± 0.0002          –0.01151 ± 0.00004         W UMa
		                     R          0.0053 ± 0.0002      –0.0595 ± 0.0002        –0.0008 ± 0.0002          –0.01097 ± 0.00004         W UMa

		                     B          0.0812 ± 0.0006      –0.1567 ± 0.0006        –0.0716 ± 0.0006          –0.0441 ± 0.0002           Algol
  MU Leo               V          0.0734 ± 0.0005      –0.1546 ± 0.0005        –0.0759 ± 0.0006          –0.0432 ± 0.0002           Algol
		                     R          0.0582 ± 0.0006      –0.1526 ± 0.0006        –0.0779 ± 0.0007          –0.0424 ± 0.0002           Algol




Figure 1. Normalized flux for LX Leo (left), V354 UMa (center), and MU Leo (right) in all three filters. The Fourier fits (continuous curves) are plotted along with
the blue, green, and red curves corresponding to B, V, and R filters, respectively. The average error in the flux for each measurement is approximately 0.003 and
0.005 for LX Leo and V354 UMa, respectively, in all filters, and about 0.0002 for MU Leo in all filters. Error bars are not shown for the sake of clarity.


function in AIJ to determine the photometric aperture radius                        3. Analysis
and the radii of the inner and outer annulus, using the method
outlined in Collins et al. (2017) (see their Appendix A). To                            The analysis of the light curves in this paper closely follows
search for comparison and check stars, we used the SIMBAD                           the procedure outlined in Gardner et al. (2015), Akiba et al.
Astronomical Database (Wenger et al. 2000) to find any stars                        (2019), and Hahs et al. (2020). We first phase-fold the time axis
in the image frame relatively close in brightness and size to our                   and ensure that the primary (deeper) eclipse always coincides
target star. We then used the cataloged B and V magnitudes of                       with phase “0.” Additionally, we calculate the normalized flux
the comparison stars to determine the corresponding magnitudes                      for each data point from the measured magnitudes obtained via
of each of our target stars (Table 1). Instrumental magnitudes                      differential photometry (Warner and Harris 2006) as:
were used for the R-filter since the R magnitude was not listed.
We also ensured that the comparison and check stars chosen                                                  I(Φ)obs = 10–0.4 × (m(Φ) – m(max))       (1)
showed no variability in each of the filters. All differential
photometry data can be retrieved from the AAVSO International                       where m(Φ) is the magnitude at a certain phase Φ and m(max)
Database (Kafka 2021). These data are also available on                             is the maximum magnitude observed for the object. We perform
request via email: gokhale@truman.edu and are available at                          Fourier fit analyses on the light curves of each object in each
http://gokhale.sites.truman.edu/asymmetries/.                                       filter similar to Wilsey and Beaky (2009). A truncated twelve-
                                                                                    term Fourier fit is given by
Hoffman et al.,       JAAVSO Volume 50, 2022                                                                                                                      5
                         12                                                        3.2. Quantifying the asymmetries in the light curves
      I(Φ)fit = a0 +     Σ (a cos(2π n Φ) + b sin(2π n Φ))   (2)
                         n=1
                               n                  n
                                                                                        We quantify the asymmetries in several different ways:
                                                                                        1. For each object, in each filter, we determine the difference
where a0, an, and bn are the Fourier coefficients of the fit, and Φ                in the normalized flux near the primary and secondary maxima
is the phase (Hoffman et al. 2009). The light curves of the three                  from the data (ΔIave ) and from the fit (ΔIfit ). Additionally, we
objects in each filter along with their Fourier fits are presented                 determine |2b1| from the Fourier coefficient b1, since this term
in Figure 1.                                                                       represents the half-amplitude of the sine wave of the Fourier
                                                                                   fit and is thus a measure of the difference between the primary
3.1. Classification                                                                and secondary maxima (Wilsey and Beaky 2009). These values
      Following Akiba et al. (2019), the Fourier coefficients                      are tabulated in Table 3. It is clear that V354 UMa has the
and the associated errors are extracted from the Fourier fits                      smallest, if any, asymmetry in its light curve, followed by MU
generated using Mathematica (Wolfram Research Co. 2019)                            Leo, for which ΔI is negative. LX Leo has the greatest amount
and are tabulated in Table 2. The values of these coefficients                     of asymmetry. Generally, ΔI is greatest in the B filter, which is
are determined by the shape of the light curves, and hence are                     consistent with the results obtained by Akiba et al. (2019).
a quantitative measure of the geometry of the eclipsing binary                          2. For each object, in each filter, we determine the O'Connell
(Ruciński 1973, 1993, 1997). In particular, the condition a4 <                     Effect ratio (OER) and the Light Curve Asymmetry (LCA,
a2 (0.125 – a2 ) implies that the system is a detached system and                  McCartney 1999) using:
is classified as an Algol-type system (see Akiba et al. 2019)
                                                                                                                  0.5
for details). If this condition is not met, and additionally if                                         ∫0.0 (I(Φ)fit – I(0.0)fit)dΦ
|a1| < 0.05, then the system is classified as a W UMa type system;                                OER = 1.0
                                                                                                       ——————————                                               (3)
                                                                                                        ∫0.5 (I(Φ)fit – I(0.0)fit)dΦ
otherwise it is classified as a β Lyrae type system.
      For the three systems under consideration, in all three filters,             and,
the results from the Fourier coefficient method are consistent


                                                                                                      √
with those from qualitative visual inspection: LX Leo and                                                   0.5
                                                                                                               (I(Φ)fit – I(1.0 – Φ)fit) 2

V354 UMa are confirmed to be W UMa type systems, while                                       LCA =        ∫0.0 —————————                 dΦ                  (7) (5)
                                                                                                                        I(Φ)fit
                                                                                                                            2
MU Leo is of the Algol-type.

Table 3. Quantifying the O’Connell Effect in terms of difference in maxima         Table 4. Quantifying the O’Connell Effect in terms of OER and LCA (see
(see section 3.2 for details).                                                     section 3.2 for the definitions of the OER and LCA).

   Target       Filter         |2b1|       ΔI (Fourier)      ΔI (Average)                    Target        Filter          OER                   LCA

		                B       0.043 ± 0.001    0.058 ± 0.007    0.053 ± 0.001          		                        B          1.08 ± 0.02          0.025 ± 0.003
 LX Leo           V       0.035 ± 0.001    0.044 ± 0.006    0.040 ± 0.001            LX Leo                  V          1.07 ± 0.02          0.020 ± 0.003
		                R       0.023 ± 0.001    0.033 ± 0.005    0.029 ± 0.001          		                        R          1.05 ± 0.02          0.014 ± 0.002

		                B      0.0020 ± 0.0002 0.000 ± 0.002      0.000 ± 0.002          		                        B          1.02 ± 0.02          0.003 ± 0.001
 V0354 UMa        V      0.0002 ± 0.0003 0.000 ± 0.002      0.000 ± 0.003            V0354 UMa               V          1.00 ± 0.03          0.002 ± 0.001
		                R      0.0016 ± 0.0003 –0.001 ± 0.002     0.000 ± 0.002          		                        R          1.02 ± 0.03          0.003 ± 0.001

		                B       0.018 ± 0.001   –0.022 ± 0.005 –0.0201 ± 0.0001          		                        B          0.98 ± 0.01          0.021 ± 0.002
 MU Leo           V       0.011 ± 0.001   –0.012 ± 0.004 –0.0123 ± 0.0001            MU Leo                  V          0.99 ± 0.01          0.013 ± 0.002
		                R       0.013 ± 0.001   –0.001 ± 0.005 –0.0120 ± 0.0001          		                        R          0.99 ± 0.01          0.015 ± 0.002




Figure 2. Difference in normalized flux between the two halves of the light curves, in the B (blue solid curve), V (green dashed), and R (red dotted) filters for
LX Leo (left plot), V354 UMa (middle), and MU Leo (right). In the absence of any asymmetry, the two curves should coincide, and the solid blue curve in the
bottom panel would be a flat line at “0.”
6                                                                                             Hoffman et al.,    JAAVSO Volume 50, 2022
The asymmetry reflected by the OER and LCA values (Table 4)              evident in the “half-phase plots” in Figure 2. In all three filters,
are consistent with the ΔI values from Table 3. We note that             the values of ΔI for V354 UMa are essentially zero implying a
V354 UMa has a very symmetric light curve, while LX Leo                  very small O’Connell effect. Similarly, the OER is quite small
has a significant asymmetry as quantified by both the OER                for V354 UMa, and the LCA is almost an order of magnitude
and LCA. Interestingly, MU Leo has a low OER in each filter,             smaller than the other two objects, implying an very symmetric
though the LCA is significant (see below). We do not observe             curve for V354 UMa.
any obvious trend in the amount of asymmetry as a function                   We do not discern any obvious correlation between the type
of filter, though generally, the OER and LCA are the greatest            of filter used and any of the asymmetry parameters we have
in the B-filter which is consistent with the results obtained            discussed, except to note that the asymmetries seem to be more
by Gardner et al. (2015), Akiba et al. (2019), and Hahs et al.           pronounced in the B-filter.
2020). Of course, a much larger sample is necessary to derive                As mentioned earlier, we are currently working on several
any reliable conclusions regarding the filter-dependance of light        EBs from the Kepler and TESS catalogs (Knote et al. 2022).
curve asymmetries in EB light curves.                                    Our goal is to extend the work presented here to hundreds of
      3. Additionally, we superpose the two halves of each of            objects, and study the time evolution of the asymmetries over
the light curves to generate “half-phase plots” to visually              several hundreds and thousands of orbital cycles. We hope to
demonstrate the asymmetries. We do this by “mirroring” the               discern patterns in the changes in the asymmetries by addressing
light curve about the phase 0.5, and in Figure 2 we plot the             questions such as: What are the timescales over which the
difference between the flux at equivalent phases in the light            asymmetries change? Are there differences in the asymmetries
curve (for example, phases Φ = 0.2 and Φ = 0.8). Consequently,           in over-contact, semi-detached, and detached systems? Why do
this plot helps to indicate the phases around which the curve is         the asymmetries in certain systems remain very steady, while in
most asymmetric, and where, for example, star spots may be               other systems the asymmetries vary significantly over relatively
located under the starpot model to explain these asymmetries. It         short timescales?, and so on.
is again clear that V354 UMa has a very symmetric light curve.
It is interesting to note that for LX Leo, in all filters, the primary   5. Acknowledgements
half of the light curve has a greater flux than the secondary half
(ΔI(Φ)fit > 0). On the other hand, in MU Leo (and to some extent             We have made extensive use of the tools available on the
in V354 UMa), the contribution to the flux varies with phase -           AAVSO website, in particular the VSP (Variable Star Plotter)
starting off with greater flux on the secondary side (ΔI(Φ)fit > 0       tool to generate star charts. In addition, we have used the
between phases 0.0–0.3, or equivalently, phases 0.7–1.0)                 SIMBAD database, operated at CDS, Strasbourg, France, and
and moving toward higher contributions from the primary                  NASA’s Astrophysics Data System. We are thankful for the
side (ΔI(Φ)fit < 0 between phases 0.3–0.5 or phases 0.5–0.7).            support provided by the Office of Student Research at Truman
Thus, even though MU Leo has an OER ≈ 1, the “half-phase”                State University and the Missouri Space Grant Consortium.
plot (Figure 2) demonstrates that there is still significant
asymmetry in its light curve. The “half-phase” plots are a               References
visual demonstration of the LCA parameter, and the example
of MU Leo demonstrates the advantages of quantifying the                 Akiba, T., Neugarten, A., Ortmann, C., and Gokhale, V. 2019,
asymmetries in different ways.                                               J. Amer. Assoc. Var. Star Obs., 47, 186.
                                                                         Collins, K. A., Kielkopf, J. F., Stassun, K. G., and Hessman, F.
4. Discussion                                                                V. 2017, Astron. J., 153, 77.
                                                                         Gardner, T., Hahs, G., and Gokhale, V. 2015, J. Amer. Assoc.
    We have quantified the asymmetries in three short period                 Var. Star Obs., 43, 186.
eclipsing binary systems: LX Leo (P = 0.235247 d), V345                  Gürol, B., Michel, R., and Gonzalez, C. 2017, Rev. Mex.
UMa (P = 0.293825 d), and MU Leo (P = 0.388442 d). Of                        Astron. Astrofis., 53, 179.
these, V354 UMa exhibits the most symmetric light curve,                 Hahs, G., Ortmann, C., and Gokhale, V. 2020, J. Amer. Assoc.
while LX Leo is the most asymmetric. We have confirmed that                  Var. Star Obs., 57, 48.
LX Leo and V354 UMa are W UMa-type systems, while MU                     Hoffman, D. I., Harrison, T. E., and McNamara, B. J. 2009,
Leo is an Algol-type system. For LX Leo, the asymmetry is                    Astron. J., 138, 466.
greatest around phase 0.25, while for MU Leo, it is greatest             Kafka, S. 2021, Observations from the AAVSO International
near phase 0.15. By visual inspection, it is clear that LX Leo               Database (https://www.aavso.org/data-download).
has a positive O’Connell Effect (primary maxima is higher than           Knote, M. K., et al. 2022, Astrophys. J., Suppl. Ser., accepted
secondary maxima) in all three filters. This bears out in the                (https://arxiv.org/abs/2206.04142).
positive values of ΔI and the fact that OER > 1. For MU Leo,             Kreiner, J. M. 2004, Acta Astron., 54, 207.
ΔI is negative in all filters, which is consistent with the shape of     McCartney, S. A. 1999, Ph.D. dissertation, University of
the light curve for MU Leo (Figure 1). Note also, for MU Leo,                Oklahoma.
the OER <  ~ 1 in all three filters, showing there is slightly less      Michel, R., Xia, Q.-Q., and Higuera, J. 2019, Res. Astron.
flux in the primary half of the light curve than the secondary               Astrophys., 19, 99 (DOI: 10.1088/1674-4527/19/7/99).
half, but overall, there is almost equal flux in the two halves.         O’Connell, D. J. K. 1951, Publ. Riverview Coll. Obs., 2, 85.
Despite this, the LCA is significant, and clear asymmetries are          Prša, A., et al. 2011, Astron. J., 141, 83.
Hoffman et al.,   JAAVSO Volume 50, 2022                                                                                 7
Ricker G. R., et al. 2015, J. Astron. Telesc. Instrum. Syst., 1,   Wilsey, N. J., and Beaky M. M. 2009, in The Society for
   014003.                                                            Astronomical Sciences 28th Annual Symposium on
Ruciński, S. M. 1973, Acta Astron., 23, 79.                           Telescope Science, Society for Astronomical Sciences,
Ruciński, S. M. 1993, Publ. Astron. Soc. Pacific, 105, 1433.          Rancho Cucamonga, CA, 107.
Ruciński, S. M. 1997, Astron. J., 113, 407.                        Wolfram Research Co. 2019, How to Fit Models with
Warner, B. D., and Harris, A. W. 2006, A Practical Guide to           Measurement Errors (https://reference.wolfram.com/
   Lightcurve Photometry and Analysis, Springer, New York.            language/howto/FitModelsWithMeasurementErrors).
Wenger, M., et al. 2000, Astron. Astrophys., Suppl. Ser., 143, 9
   (DOI: 10.1051/aas:2000332).

```
