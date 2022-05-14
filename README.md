
# BIM_processing

This program is a collaboration with Fengyan Zhang(zmocheng@gmail.com), I am mainly responsible for the Convex Hull algorithm, Fengyan is responsible for the algorithm of splitting into individual OBJs and transforming into Polyhedra. Thanks Fengyan!

This program is to implement a automatic BIM(IFC model) to Geo(CityJSON model) conversion for one building. It is implemented in C++ with CGAL library, and Niels Lohmann’s JSON for Modern C++ library.

The workflow of the program is illstrated below:
![workflow](https://user-images.githubusercontent.com/75926656/168439483-1f263293-7644-4bc6-872a-c915ae80a4c7.png)


the original Ifc file looks like:

<img width="400" alt="kit file" src="https://user-images.githubusercontent.com/75926656/168441435-8f789f56-4332-456d-94c9-62328e9992c5.png"> <img width="400" alt="kit file2" src="https://user-images.githubusercontent.com/75926656/168441442-eb62b196-1d33-4c77-8d4f-4fe926246c84.png">

We use IfcConvert to convert IFC model to OBJ model:

<table><tr><td bgcolor=gray>IfcConvert --include+=entities IfcSlab IfcDoor IfcWindow IfcWall --weld-vertices --orient-shells --validate KIT.ifc KIT.obj</td></tr></table>

(KIT.ifc is the self-defined filename)

The outer shell and the inner shell looks like:

<img width="400" alt="shell" src="https://user-images.githubusercontent.com/75926656/168441452-503c65ea-e539-4365-8341-2e74e6bdc49d.png"> <img width="400" alt="inner_shell" src="https://user-images.githubusercontent.com/75926656/168441456-231b4856-459b-4676-884f-161a36876e5a.png">

tips: this program is more suitable for buildings with only one floor.
