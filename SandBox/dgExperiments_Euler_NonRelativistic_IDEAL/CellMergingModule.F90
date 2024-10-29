MODULE CellMergingModule

    USE KindModule, ONLY: &
        DP
    USE GeometryFieldsModule, ONLY: &
        uGF, iGF_h_2, iGF_h_3
    USE MeshModule, ONLY: &
        MeshX
    USE LagrangePolynomialsModule, ONLY: &
        LagrangeP

    PUBLIC :: Initialize_CellMerging
    PUBLIC :: Finalize_CellMerging
    PUBLIC :: Determine_MergedCells
    PUBLIC :: MergetoCoarseGrid
    PUBLIC :: RestricttoFineGrid

CONTAINS

    SUBROUTINE Initialize_CellMerging

    END SUBROUTINE Initialize_CellMerging

    SUBROUTINE Finalize_CellMerging

    END SUBROUTINE Finalize_CellMerging

    SUBROUTINE Determine_MergedCells

    END SUBROUTINE Determine_MergedCells

    SUBROUTINE MergetoCoarseGrid

    END SUBROUTINE MergetoCoarseGrid

    SUBROUTINE RestricttoFineGrid

    END SUBROUTINE RestricttoFineGrid

END MODULE CellMergingModule