using PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using Cbc
using Dates
const PSY = PowerSystems
using CPLEX
using PowerGraphics
using JuMP
plotlyjs()

include("get_data.jl")

free_solver = optimizer_with_attributes(Cbc.Optimizer, "ratioGap" => 1e-3)
solver = optimizer_with_attributes(CPLEX.Optimizer)

# Custom Generator
mutable struct ThermalMustRun <: PSY.ThermalGen
    name::String
    available::Bool
    status::Bool
    bus::PSY.Bus
    activepower::Float64
    reactivepower::Float64
    rating::Float64
    primemover::PrimeMovers.PrimeMover
    fuel::ThermalFuels.ThermalFuel
    activepowerlimits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    reactivepowerlimits::NamedTuple{(:min, :max),Tuple{Float64,Float64}}
    op_cost::TwoPartCost
    basepower::Float64
    internal::PSY.InfrastructureSystemsInternal
end

PSY.get_name(d::ThermalMustRun) = d.name
PSY.get_available(d::ThermalMustRun) = d.available
PSY.get_bus(d::ThermalMustRun) = d.bus
PSY.get_status(d::ThermalMustRun) = d.status
PSY.get_activepower(d::ThermalMustRun) = d.activepower
PSY.get_reactivepower(d::ThermalMustRun) = d.reactivepower
PSY.get_rating(d::ThermalMustRun) = d.rating
PSY.get_primemover(d::ThermalMustRun) = d.primemover
PSY.get_fuel(d::ThermalMustRun) = d.fuel
PSY.get_activepowerlimits(d::ThermalMustRun) = d.activepowerlimits
PSY.get_reactivepowerlimits(d::ThermalMustRun) = d.reactivepowerlimits
PSY.get_op_cost(d::ThermalMustRun) = d.op_cost
PSY.get_basepower(d::ThermalMustRun) = d.basepower

sys_DA = System("RTS_1hr_sys_agc.json")
bus = get_component(Bus, sys_DA, "Bacon")

mythermal = ThermalMustRun(
    "Pickering",
    true,
    true,
    bus,
    4.0,
    0.0,
    5.08,
    PrimeMovers.ST,
    ThermalFuels.NUCLEAR,
    (min = 3.9, max = 4.1),
    (min = -3.0, max = 3.0),
    TwoPartCost(0.0, 0.0),
    1.0,
    PSY.InfrastructureSystemsInternal()
)

add_component!(sys_DA, mythermal)

#### Use my custom generator with an available formulation.
# Caveats: The use of available formulations require that certain getter functions are included
# and return valid quantitites.

uc_problem = UnitCommitmentProblem(sys_DA; initial_time = DateTime("2020-09-01"), horizon = 48)

MustRunModel = DeviceModel(ThermalMustRun, ThermalDispatch)

construct_device!(uc_problem, :MustRun, MustRunModel)

results_uc = solve!(uc_problem, optimizer = solver)
fuel_plot(results_uc, sys_DA)
stack_plot(results_uc)

#### Developing my custom model for my custom type #####

struct MustRunFormulation <: PSI.AbstractThermalFormulation end

### Customization Option 1 (Most Complicated, try not to do this):
# Fully Custom variable and constraint addition

function PSI.construct_device!(op_problem::OperationsProblem, label::Symbol, model::DeviceModel{ThermalMustRun, MustRunFormulation})

    time_steps = PSI.model_time_steps(op_problem.psi_container)
    devices_of_my_type = PSY.get_components(ThermalMustRun, op_problem.sys)
    jump_model = PSI.get_jump_model(op_problem)

    # If these steps are skipped, the variable won't be added to the psi_container
    variable_container = PSI.container_spec(jump_model,
        (PSY.get_name(d) for d in devices_of_my_type),
        time_steps)
    PSI.assign_variable!(op_problem.psi_container, PSI.variable_name(ACTIVE_POWER, ThermalMustRun), variable_container)

    # Add variable to container and to the expression of the nodal balance
    for t in time_steps, d in devices_of_my_type
        variable_container[get_name(d),t] = JuMP.@variable(jump_model, lower_bound = 0.0)
        bus_number = PSY.get_number(PSY.get_bus(d))
        PSI.add_to_expression!(
                PSI.get_expression(op_problem.psi_container, :nodal_balance_active),
                bus_number,
                t,
                variable_container[get_name(d),t],
                1.0,
            )
    end

    # Create the constraint
    constraint_container = PSI.JuMPConstraintArray(undef, (PSY.get_name(d) for d in devices_of_my_type), time_steps)
    PSI.assign_constraint!(op_problem.psi_container, "fix_MustRun", constraint_container)
    for t in time_steps, d in devices_of_my_type
       name = get_name(d)
       constraint_container[name, t] = JuMP.@constraint(jump_model,
                                        variable_container[name, t] == get_activepower(d))
    end

    return
end

uc_problem = UnitCommitmentProblem(sys_DA; initial_time = DateTime("2020-09-01"), horizon = 48)

MustRunModel = DeviceModel(ThermalMustRun, MustRunFormulation)

construct_device!(uc_problem, :MustRun, MustRunModel)

results_uc = solve!(uc_problem, optimizer = solver)
stack_plot(results_uc)

### Customization Option 2 (Use modeling low level provided functions):
# Fully Custom variable and constraint addition Currently being improved

function PSI.construct_device!(op_problem::OperationsProblem, label::Symbol, model::DeviceModel{ThermalMustRun, MustRunFormulation})
    time_steps = PSI.model_time_steps(op_problem.psi_container)
    devices_of_my_type = PSY.get_components(ThermalMustRun, op_problem.sys)
    jump_model = PSI.get_jump_model(op_problem)

    PSI.add_variable(
        op_problem.psi_container,
        devices_of_my_type,
        PSI.variable_name(ACTIVE_POWER, ThermalMustRun),
        false,
        :nodal_balance_active;
        ub_value = d -> PSY.get_activepowerlimits(d).max,
        lb_value = d -> PSY.get_activepowerlimits(d).min,
        init_value = nothing,
    )

    PSI.add_variable(
        op_problem.psi_container,
        devices_of_my_type,
        PSI.variable_name(ON, ThermalMustRun),
        false,
        :nodal_balance_active;
        ub_value = d -> 1.0,
        lb_value = d -> 1.0,
        init_value = 1.0,
    )

    PSI.activepower_constraints!(op_problem.psi_container,
                                 devices_of_my_type,
                                 model, op_problem.template.transmission, nothing)

    PSI.cost_function(op_problem.psi_container, devices_of_my_type, ThermalDispatch, op_problem.template.transmission, nothing)

    return
end

uc_problem = UnitCommitmentProblem(sys_DA; initial_time = DateTime("2020-09-01"), horizon = 48)

MustRunModel = DeviceModel(ThermalMustRun, MustRunFormulation)

construct_device!(uc_problem, :MustRun, MustRunModel)

results_uc = solve!(uc_problem, optimizer = solver)
stack_plot(results_uc)


### Making a custom opertaions problem #####

# 1. You need your custom OperationsProblem

struct MustRunOperationProblem <: PSI.AbstractOperationsProblem end

# 2. Your template

MustRunModel = DeviceModel(ThermalMustRun, MustRunFormulation)
my_template = template_unit_commitment()
PSI.set_model!(my_template, :MustRun, MustRunModel)
PSI.set_model!(my_template, DCPPowerModel)

problem = OperationsProblem(MustRunOperationProblem, my_template, sys_DA; initial_time = DateTime("2020-09-01"), horizon = 48)
