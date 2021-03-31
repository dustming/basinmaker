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
  