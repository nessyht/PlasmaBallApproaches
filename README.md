<h1>Simulation Approaches for a Virtual Plasma Globe</h1>

An advanced practical project 

Exploration of some simulation approaches for a Virtual Plasma Globe as a precursor for a thesis.

The task of relatively heavy iterative simulation suited a parallel approach. The software requires OpenCL 1.2, it can be reconfigured to run as a multi-threaded application on CPU if no discrete device is available.

  
<h2>1 Introduction</h2>
<img width="30%" align="right" alt="Shows an a photograph of a Plasma ball." src="https://user-images.githubusercontent.com/44576195/187530672-32a3c359-6b72-4bdf-b1f3-c2db5a2888ae.jpg"><p>Plasma globes are commercial toys that create
glowing filaments in a hollow glass ball. These flalments "dance" in a light display caused by electric
arc discharge in a gas at atmospheric pressure. Visually, the filaments seem stable for a short time,
as though fixed to end-points at the central electrode and at the outer edge of the glass. These
endpoints drift along the outer surface for a few
centimetres before the filament disappears. Many
of the filaments branch out to reach the outer glass
at multiple points per filament. Some numerical
models for electric discharge exist. The aim of
this project, is to explore the application and limitations of the dielectric breakdown model to simulate the filaments within a plasma globe. This
will involve expanding the model to include multiple filaments, to test how the filaments interact
with each other and to investigate their movement.
Ultimately the aim is to test the effcacy of the
Dielectric Breakdown model as the core tool for
simulating and visualising a virtual plasma ball as
a proof of concept for future work.</p>
  
Full report can be found as pdf: 
<a href="https://github.com/nessyht/Plasma_globe/blob/67c3aed0963ebbf113d05ea2609aab0af67baa1f/LaTeX/main.pdf">SimulationApproachesForAVirtualPlasmaGlobe</a>
</br>
Along with the accompanying presentation as pptx:
<a href="https://github.com/nessyht/Plasma_globe/blob/8336d9783d04dfad3f5cca3d576c8345927dfe11/LaTeX/Simulation%20Approaches%20for%20a%20Virtual%20Plasma%20Globe.pptx">SimulationApproachesForAVirtualPlasmaGlobe</a>
</br>

<h2>Snapshot - Results</h2>

<img width="23%" align="center" alt="Shows an a photograph of a Plasma ball." src="https://user-images.githubusercontent.com/44576195/187530799-64475515-b6a9-4bc9-9e0d-725e3fb82e62.png"> <img width="23%" align="center" alt="Shows an a photograph of a Plasma ball." src="https://user-images.githubusercontent.com/44576195/187530800-d76befb2-d6ce-458a-8c4f-ab290d123d37.png"> <img width="23%" align="center" alt="Shows an a photograph of a Plasma ball." src="https://user-images.githubusercontent.com/44576195/187530804-902da82b-a5a1-468b-9f3b-9c5f0377bc05.png"> <img width="23%" align="center" alt="Shows an a photograph of a Plasma ball." src="https://user-images.githubusercontent.com/44576195/187530806-468e5eaa-70f5-42e3-a116-cc0b8c8c788f.png">
<caption>Snapshots of simulated filaments in a virtual plasma ball of radius 128-cells<caption>

Results are generated as a series of .vtk files, these can be opened using Paraview or similar programs
