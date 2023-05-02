## INTRODUCTION: 
Optimal Design of Compact One-Stage Spur Gear Reducer using Genetic Algorithm. Code was written as part of a class project for ME 4244: Machine Design II taught by Dr. Md. Wahah in Spring 2023 in Louisiana State University. 

Written in MATLAB. Tested in MATLAB 2022a. For easier debugging and future upgrades, I have split important functions into their own .m files

The system finds the optimal values of number of teeth for pinion gear (driving gear), pitch module and and face width. All units in Imperial units.Objective function defined as minimization of total system weight where total weight is calculated as sum of two gears weight.

## HOW TO USE:
+ Use the GEAR_DESIGN_MAIN_NOTEBOOK.mlx notebook. It is self contained and will guide you through each step 

## NOTES
+ Genetic Algorithm code based on https://github.com/Mechazo11/Genetic_Algorithm_2D_Rosenbrock

+ Fatigue analysis based on Gear-Tooth Bending Strength and Gear-Tooth Surface Fatigue Analysis Methodologies described in Chapter 15 "Spur Gears", Fundamental of Machine Component Design, 7E by Robert C. Juvinall and Kurt M. Marshek.

+ Objective function defined as a weight minimization problem.

+ Gears designed for 10^6 cycles (i.e infinite life)


## ADDITIONAL MATERIALS 
+ A pdf copy of that Chapter 15
+ A pdf copy of the report that used this code for the class project
+ Some of the papers  

## DISCLAIMER
THIS SOFTWARE WAS DEVELOPED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE USED FOR COMMERCIAL PURPOSES. THE SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE (ChatGPT wrote it :)).
