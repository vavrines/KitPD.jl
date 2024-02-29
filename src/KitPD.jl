module KitPD

using KitBase
using KitBase.Parameters: @with_kw
import Peridynamics as PD

export KP, PD
export PDMater

const KP = KitPD

include("type.jl")

end
