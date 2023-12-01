//==========================================================================
//  Implementation of KLM Calorimeter
//--------------------------------------------------------------------------
//  Author: Sebouh Paul (UCR)
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml_det_t     detElem = handle;
  std::string   detName = detElem.nameStr();
  int           detID   = detElem.id();

  xml_dim_t  dim        = detElem.dimensions();
  double     radius     = dim.x(); // Size in radius
  double     length     = dim.z(); // Size along z-axis

  xml_dim_t  pos        = detElem.position(); // Position in global coordinates

  Material   air        = desc.material("Air");
  
  // Getting beampipe hole dimensions
  const xml::Component &insert_hole_xml = detElem.child(
    _Unicode(insert_hole)
  );
  const double hole_height = dd4hep::getAttrOrDefault<double>(
    insert_hole_xml,
    _Unicode(hole_height),
    60*cm
  );
  const double hole_width = dd4hep::getAttrOrDefault<double>(
    insert_hole_xml,
    _Unicode(hole_width),
    60*cm
  );

  const double hole_offset_x = dd4hep::getAttrOrDefault<double>(
    insert_hole_xml,
    _Unicode(hole_offset_x),
    -10*cm
  );

    
  // Defining envelope
  Tube envelope(0, radius, length/2);

  /*
    Hole initial x-position in local coordinate system
    x = 0 in local coordinate system is x = HcalEndcapPInsert_xposition in global coordinate system
    Hole starts at x = 0 cm in global system
  */
  //const double hole_x_1 = 0 - pos.x();

  Box insert_hole(hole_width/2., hole_height/2., length / 2.);
  SubtractionSolid envelope_with_hole(
	envelope,
        insert_hole,
        Position(
		 hole_offset_x,
          0.,
          0
        )
      );
  //commented out because it was causing print_materials to crash
  //envelope = envelope_with_hole;
    

  // Defining envelope volume
  Volume envelopeVol(detName, envelope, air); 
  // Setting envelope attributes
  envelopeVol.setAttributes(
    desc,
    detElem.regionStr(),
    detElem.limitsStr(),
    detElem.visStr()
  );

  PlacedVolume pv;

  
  double z_distance_traversed = 0.;
  int layer_num = 1;
  
  // Looping through all the different layer sections
  for(xml_coll_t c(detElem,_U(layer)); c; ++c)
  {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();  
    double layer_thickness = x_layer.thickness();
    
    // Looping through the number of repeated layers in each section
    for(int i = 0; i < repeat; i++)
    {
      std::string layer_name = detName + _toString(layer_num, "_layer%d");

      Tube layer(0, radius, layer_thickness / 2.);

      // Removing insert shape from each layer

      IntersectionSolid layer_with_hole(
        layer,
        envelope
      );
      Volume layer_vol(layer_name, layer_with_hole, air);

      

      int slice_num = 1;
      double slice_z = -layer_thickness / 2.; // Keeps track of slices' z locations in each layer

      // Looping over each layer's slices
      for(xml_coll_t l(x_layer,_U(slice)); l; ++l) 
      {
        xml_comp_t x_slice = l;
        double slice_thickness = x_slice.thickness();
        std::string slice_name = layer_name + _toString(slice_num, "slice%d");
        Material slice_mat = desc.material(x_slice.materialStr());		
        slice_z += slice_thickness/2.; // Going to slice halfway point

        // Each slice within a layer has the same hole radius
        Tube slice(0, radius, slice_thickness/2.);
        
        IntersectionSolid slice_with_hole(
          slice,
          envelope
        );
        Volume slice_vol (slice_name, slice_with_hole, slice_mat);
        
        // Setting appropriate slices as sensitive
        if(x_slice.isSensitive())
        {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }
        
        // Setting slice attributes
        slice_vol.setAttributes(
          desc,
          x_slice.regionStr(),
          x_slice.limitsStr(),
          x_slice.visStr()
        );

        // Placing slice within layer
        pv = layer_vol.placeVolume(
          slice_vol,
          Transform3D(
            RotationZYX(0, 0, 0),
            Position(
              0.,
              0.,
              slice_z
            )
          )
        );
        pv.addPhysVolID("slice", slice_num);
        slice_z += slice_thickness/2.;
        z_distance_traversed += slice_thickness;
        ++slice_num;
      }

      // Setting layer attributes
      layer_vol.setAttributes(
        desc,
        x_layer.regionStr(),
        x_layer.limitsStr(), 
        x_layer.visStr()
      );
      // Placing each layer inside the envelope volume
      // -length/2. is front of detector in global coordinate system
      // + (z_distance_traversed - layer_thickness) goes to the front of each layer
      // + layer_thickness/2. places layer in correct spot
      // Example: After placement of slices in first layer, z_distance_traversed = layer_thickness
      //          Subtracting layer_thickness goes back to the front of the first slice (Now, z = -length/2)
      //          Adding layer_thickness/2. goes to half the first layer thickness (proper place to put layer)
      //          Each loop over repeat will increases z_distance_traversed by layer_thickness
      pv = envelopeVol.placeVolume(
        layer_vol,
        Transform3D(
          RotationZYX(0, 0, 0),
          Position(
            0.,
            0.,
            -length/2.  + (z_distance_traversed - layer_thickness) + layer_thickness/2.
          )
        )
      );
      pv.addPhysVolID("layer", layer_num);
      layer_num++;
      //hole_x_tracker -= 0.0; // The hole shifts along -x by 0 cm every layer
    }
  }

  DetElement   det(detName, detID);  
  Volume motherVol = desc.pickMotherVolume(det);
  
  // Placing insert in world volume
  auto tr = Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.));
  PlacedVolume phv = motherVol.placeVolume(envelopeVol, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(ForwardHcalKLM, createDetector)
