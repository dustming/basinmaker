========
Overview
========

Before using the BasinMaker or routing product, the differences between routing structure without considering lake and routing structure with lakes defined by BasinMaker and routing product will be introduce in this section. 

A routing structure without considering lakes is showed in Figure A. Catchments in a routing structure without considering lakes are only determined by river reaches (Figure A). When a semi hydrological model build with this routing structure, the lakes inflow and outflow cannot be simulated, and thus the impact of lake on the routing process modeling such as flow attenuation cannot be correctly modeled. The reason is that the streamflow is only explicitly simulated at the outlet of each catchment in semi-distributed hydrological models. It is not explicitly simulated at the mid of the river segment or any point inside each catchment. 
An example of routing structure with lakes represented by the BasinMaker or routing product is shown in Figure B. Lakes are divided into two categories: (1) connected lakes (CL), which indicates lakes outlets are explicitly connected to a downstream non-zero length river channel/lake in the routing product; and (2) non-connected lakes (NCL), which denotes lakes is not explicitly connected to the downstream routing network. The connection is more implicit (Figure B). 
Both connected lakes (CLs) and NCLs within a watershed are considered to be contributing areas of the watershed. As such, both CLs and NCLs will drain to the outlet of the watershed. Both CLs and NCLs are represented by a lake catchment (Figure B). A lake catchment is defined by the following rules:1) The extent of the lake catchment will fully cover the lake; 2) the outlet of the lake catchment is the same as the outlet of the lake; 3) each lake’s inlets are treated as a catchment outlet. In this way, both inflow and outflow of each lake can be explicitly simulated by hydrologic routing models.

The only difference between CLs and NCLs is that CLs always drain into an explicitly represented river channel that is connected to the lake outlet while NCLs do not. NCLs exist because for smaller catchments, flow accumulation threshold settings can sometimes suppress the creation of a river channel at the lake outlet. As such, NCLs should drain directly, via a zero length flow path, into the next downstream catchment. Users need to ensure their hydrologic routing model accomplishes this. Specific hydrologic routing logic is as follows for NCLs: NCL catchment outflows need to be delivered to the next downstream river channel and if that river channel has a zero length (only possible if downstream catchment is also a lake catchment), that water must be delivered directly to the lake in this downstream catchment. 

.. image:: https://github.com/dustming/RoutingTool/wiki/Figures/Figure1.png
  :width: 900
  :alt: Alternative text


========  
Routing products from BasinMaker 
========

 - The North American Lake-River Routing Product v2.1 developed using BasinMaker is available at `here <http://hydrology.uwaterloo.ca/basinmaker/index.html>`_.
 - BasinMaker generated Ontario Lake-River Routing Product v1 is available at `here <https://lake-river-routing-products-uwaterloo.hub.arcgis.com>`_. .

========  
BasinMaker on Google Colab
========
 
A post-processing example via google colab can be found at here `here <https://colab.research.google.com/drive/14OC8l4ZeabOGGi0bL0ZFK1QzTOY8M9yM?usp=sharing>`_. The google colab is an online python notebook dose not require installation. This example will show you how to discretize, simplify, and revise the provided routing product for your purposes. 


========  
Authors
========
  
BasinMaker and the associated river and lake routing product was developed by the hydrology research group at the University of Waterloo. Primary Contributors are Ming Han, Hongren Shen, Bryan A. Tolson, James R. Craig and Juliane Mai. Secondary contributors are Simon Lin, Nandita B. Basu and Frezer Awol. 
  
========
Citation
========  
    
- Han, M., H. Shen, B. A. Tolson, J. R. Craig, J. Mai, S. Lin, N. B. Basu, F. Awol. 2021. BasinMaker 3.0: a GIS toolbox for distributed watershed delineation of complex lake-river routing networks. Submitted in June 2022 to Environmental Modelling and Software.
  
and/or cite the zenodo dataset directly:
  
- Ming Han, Hongren Shen, Bryan A. Tolson, & Robert A. Metcalfe. (2022). Ontario Lake-River Routing Product version 1.0 (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6536085
   
- Ming Han, Hongren Shen, Bryan A. Tolson, James R. Craig, Juliane Mai, Simon Lin, Nandita Basu, & Frezer Awol. (2020). North American Lake-River Routing Product v 2.1, derived by BasinMaker GIS Toolbox (v2.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4728185
    
========  
License
========  

BasinMaker is open-source under the `Artistic License 2.0 <https://opensource.org/licenses/Artistic-2.0>`_. This sofware is freely distributed ’as is’ without warranties or conditions of any kind, either express or implied, including, without limitation, any warranties or conditions of title, non-infringement, merchantability, or fitness for a particular purpose.
