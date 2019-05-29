/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myMovingWallSlipFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myMovingWallSlipFvPatchVectorField::
myMovingWallSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::myMovingWallSlipFvPatchVectorField::
myMovingWallSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::myMovingWallSlipFvPatchVectorField::
myMovingWallSlipFvPatchVectorField
(
    const myMovingWallSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::myMovingWallSlipFvPatchVectorField::
myMovingWallSlipFvPatchVectorField
(
    const myMovingWallSlipFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::myMovingWallSlipFvPatchVectorField::
myMovingWallSlipFvPatchVectorField
(
    const myMovingWallSlipFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myMovingWallSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    if (mesh.moving())
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        const scalar deltaT = mesh.time().deltaTValue();

        const vectorField Up((pp.faceCentres() - oldFc)/deltaT);  // lg: field of wall velocity

        const volVectorField& U =
            static_cast<const volVectorField&>(internalField()); // lg: field of given vector in U input file
        scalarField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );
        
        

        const vectorField n(p.nf()); // lg: assuming this is unit surface normal field of patch
        const scalarField& magSf = p.magSf();
        tmp<scalarField> Un = phip/(magSf + VSMALL);
        
     const vectorField nHat(this->patch().nf());

        vectorField::operator=(nHat*(nHat &  (Up + n*(Un - (n & Up))) )
		+ transform(I - sqr(nHat), this->patchInternalField()) );

    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::myMovingWallSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        myMovingWallSlipFvPatchVectorField
    );
}

// ************************************************************************* //
