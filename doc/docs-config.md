src_dir: ../Modules
       ../W
       ../Chi_RPA
output_dir: ./build
exclude: Optabs.f
       ../W/PointR.f
       ../W/gjel.f
       ../W/Sloss.f
       ../W/paths.f
       ../W/f77_to_f90.f90
       ../Chi_RPA/Optabs.f
       ../Chi_RPA/Optabs.f90
       ../Chi_RPA/paths.f
       ../Chi_RPA/paths.f90
       ../Chi_RPA/paths_old.f90
       ../Chi_RPA/paths__genmod.f90
       ../Chi_RPA/pathsb__genmod.f90
       ../Chi_RPA/pathslegacy__genmod.f90
       ../Chi_RPA/PointR.f
       ../Chi_RPA/PointR_old.f
       ../Chi_RPA/PointR_old.f90
       ../Chi_RPA/gjel.f
project_github: https://github.com/Nevensky/2DMaterialsOpticalProps
summary: Fortran package for computing optical properties of quasi-2D materials.
author: Vito Despoja, Neven Golenić
email: neven.golenic@gmail.com
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: true 
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: GNU General Public License
extra_filetypes: sh #

This is the project documentation for a set of packages which compute various optical propeties of quasi-2D materials jointly developed and maintained by Vito Despoja and Neven Golenić.


@note
The code is divided in four separate programs:

- Optabs.f90 - Calculates the RPA current-current response function.
- Sloss.f90 - Calculates the GW self-energy and the screened Coulomb interaction.
- BSE.f90 - Solves the Bethe-Salpeter equation.
- Photon.f90 - Solves the Dyson equation for the photon propagator given a response function.

@endnote


## Theoretical innuendo
\begin{equation*}
\chi = 1-\chi_0 V
\end{equation*}

