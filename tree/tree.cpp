#include "tree.h"

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <math.h>

#include "../misc/utils.h"
#include "../misc/model.h"
#include <iso646.h>

using std::cout; using std::endl;
using std::stringstream; using std::string;
using std::copy; using std::fill;
using std::array;



Tree::Tree(BodyManager* bm, const double * min_bounds, const double * max_bounds, double theta){
    /* parameter for opening criterion */
    m_theta = theta;
    bodyManager = bm;
    /* create root cell */
    root = new Cell();
    copy(min_bounds, &min_bounds[3], root->min_bounds);
    copy(max_bounds, &max_bounds[3], root->max_bounds);
    root->m = 0;
    fill(root->subcells, &root->subcells[8], nullptr);
}

Tree::~Tree(){
    delete_subtree(root);
}

void Tree::delete_subtree(Cell * cell){
    /* delete subcells first */
    for(int i = 0; i < 8; i++){
        if(cell->subcells[i] != nullptr){
            delete_subtree(cell->subcells[i]);
        }
    }
    /* delete cell */
    delete cell;
}

void Tree::insert_body(int bodyIndex) {
    insert_body(root, bodyIndex);
}

void Tree::insert_body(Cell* cell, int bodyIndex)
{
    auto bodies = bodyManager->localBodies;
    /* if cell does not contain a body and is a leaf */
    if (cell->body == -1 and cell->subcells[0] == nullptr) 
    {
        cell->m = bodies.mass[bodyIndex];
    	
        copy(bodies.position.begin() + bodyIndex * 4, bodies.position.begin() + bodyIndex * 4 + 2, cell->rm);
        cell->body = bodyIndex;
    }
    else 
    {
        /* if cell already contains a body, but is not split */
        if (cell->subcells[0] == nullptr) {
            /* split the cell  */
            for (int i = 0; i < 8; i++) {
                /* create subcell */
                Cell* subcell = new Cell();

                subcell->m = 0;

                /* find min and max bounds of cell */
                for (int c = 0; c < 3; c++) {
                    double half_side = (cell->max_bounds[c] - cell->min_bounds[c]) / 2;
                    int shift = ((i >> c) & 1);
                    subcell->min_bounds[c] = cell->min_bounds[c] + shift * half_side;
                    subcell->max_bounds[c] = cell->max_bounds[c] - !shift * half_side;
                }

                /* set as subcell */
                cell->subcells[i] = subcell;

                // insert the old body into new cell  
                if (cell->body != -1 and coord_in_cell(cell->subcells[i], cell->rm)) 
                {
                    insert_body(cell->subcells[i], cell->body);
                    cell->body = -1;
                }
            }
        }

        /* insert new body into subcell */
        for (int i = 0; i < 8; i++) 
        {
            if (cell->subcells[i] != nullptr and coord_in_cell(cell->subcells[i], bodies.position.data() + bodyIndex * 4 * sizeof(double))) 
            {
                insert_body(cell->subcells[i], bodyIndex);
                break;
            }
        }

        /* update mass and center of mass for cell */
        for (int c = 0; c < 3; c++) 
        {
            cell->rm[c] = (cell->m * cell->rm[c] + bodies.mass[bodyIndex] * bodies.position[bodyIndex * 4 + c]) / (cell->m + bodies.mass[bodyIndex]);
        }
        cell->m += bodies.mass[bodyIndex];
    }
}
    
void Tree::insert_emptycell(const double * min_bounds, const double * max_bounds){
   insert_emptycell(root, min_bounds, max_bounds); 
}

void Tree::insert_emptycell(Cell * cell, const double * min_bounds, const double * max_bounds){
    /* loop through subcells */
    for(int i = 0; i < 8; i++){
        /* if bounds in subcell, recursively insert */
        if(cell->subcells[i] != nullptr and bounds_in_cell(cell->subcells[i], min_bounds, max_bounds)){
            insert_emptycell(cell->subcells[i], min_bounds, max_bounds); 
            return;
        }
        /* if bounds does not fit in subcell, create a new one */
        else if(cell->subcells[i] == nullptr){
            Cell * subcell = new Cell();
            copy(min_bounds, &min_bounds[3], subcell->min_bounds);
            copy(max_bounds, &max_bounds[3], subcell->max_bounds);
            subcell->m = 0;
            fill(subcell->subcells, &subcell->subcells[8], nullptr);
            cell->subcells[i] = subcell;
            return;
        }
    }
    throw "ERROR: Could not fit the empty cell in any cell in the tree";
}

void Tree::insert_cell(Cell * cell_to_insert){
    return insert_cell(root, cell_to_insert);
}

void Tree::insert_cell(Cell * cell, Cell * cell_to_insert){

    /* update mass and center of mass */
    if(cell->m != 0 or cell_to_insert->m != 0){
        for(int c = 0; c < 3; c++){
            cell->rm[c] = (cell->m * cell->rm[c] + cell_to_insert->m * cell_to_insert->rm[c]) / (cell->m + cell_to_insert->m);
        }
        cell->m += cell_to_insert->m;
    }

    /* insert cell in subcell */
    for(int i = 0; i < 8; i++){
        /* if bounds fit in cell */
        if(cell->subcells[i] != nullptr and bounds_in_cell(cell->subcells[i], 
                                                cell_to_insert->min_bounds, cell_to_insert->max_bounds)){
            /* recursively insert */
            insert_cell(cell->subcells[i], cell_to_insert); 
            return;
        }
        /* if bounds does not fit in any subcell add as subcell */
        else if(cell->subcells[i] == nullptr){
            cell->subcells[i] = cell_to_insert;
            return;
        }
    }
    throw "ERROR: Could not fit the cell in any cell in the tree";
}

void Tree::insert_cell(const double * min_bounds, const double * max_bounds, double m, const double * rm){
    insert_cell(root, min_bounds, max_bounds, m, rm);
}
    
void Tree::insert_cell(Cell * cell,
                    const double * min_bounds, const double * max_bounds,
                    double m, const double * rm){
    /* if a cell with same bounds already exists */
    if(same_cell(cell, min_bounds, max_bounds)){
        /* update mass and center of mass */
        if(cell->m != 0 or m != 0){
            for(int c = 0; c < 3; c++){
                cell->rm[c] = (cell->m * cell->rm[c] + m * rm[c]) / (cell->m + m);
            }
            cell->m += m;
        }
        return;
    }
 
    /* otherwise, insert cell in subcell */
    for(int i = 0; i < 8; i++){
        /* if bounds fit in cell */
        if(cell->subcells[i] != nullptr and bounds_in_cell(cell->subcells[i], min_bounds, max_bounds)){
            /* recursively insert */
            insert_cell(cell->subcells[i], min_bounds, max_bounds, m, rm); 
            return;
        }
        /* if bounds does not fit in subcell, create a new one */
        else if(cell->subcells[i] == nullptr){
            Cell * subcell = new Cell();
            copy(min_bounds, &min_bounds[3], subcell->min_bounds);
            copy(max_bounds, &max_bounds[3], subcell->max_bounds);
            copy(rm, &rm[3], subcell->rm);
            subcell->m = m;
            fill(subcell->subcells, &subcell->subcells[8], nullptr);
            subcell->inserted=true;
            cell->subcells[i] = subcell;
            return;
        }

    }
    
    throw "ERROR: Could not fit the cell in any cell in the tree";
}
    
string Tree::to_string(bool fulltree) const {
    return subtree_str(root, fulltree);
}

string Tree::subtree_str(Cell * cell, bool fulltree) const{
    /* stream to stringstream */
    std::stringstream sstream;
        
    
    int count = 0;
    if(!fulltree){
        for(count = 0; count < 8; count++){
            if(cell->subcells[count] == nullptr){
                break;
            }
        }
    }

    /* if not leaf */
    if(cell->subcells[0] != nullptr and (fulltree or count <= 2)){
        /* call function for each subcell */
        for(int i = 0; i < 8; i++){
            if(cell->subcells[i] != nullptr){
                sstream << subtree_str(cell->subcells[i], fulltree);
            }
        }
    }
    /* if leaf */
    else{
        /* output center of mass */
        //if(cell->m != 0){
        //    for(int c = 0; c < 3; c++){
        //        sstream << cell->rm[c] << " ";
        //    }
        //    sstream << std::endl;
        //}
        /* output cell boundaries */
        for(int c = 0; c < 3; c++){
            sstream << cell->min_bounds[c] << " " << cell->max_bounds[c] << " ";
        }
        sstream << std::endl;
        /* tell if cell was imported from other process */
        //if(cell->inserted){
        //     sstream << "1\n";
        //}else{
        //     sstream << "0\n";
        //}
    }

    return sstream.str();
}

bool Tree::opening_criterion(const Cell * cell, const double * min_bounds, const double * max_bounds){
    /* minimum squared distance between cell and domain */
    double d2 = dist_aabb(cell->min_bounds, cell->max_bounds, min_bounds, max_bounds);

    /* squared volume of cell */
    double v2 = pow(cell_volume(cell),0.333*2);

    /* opening criterion rule */
    if(m_theta * m_theta * d2 < v2){
        return true;
    }
    else{
        return false;
    }
}

void Tree::cells_to_send(const double * min_bounds, const double * max_bounds,
                    int min_depth, vector<Cell*> & cells){
    cells_to_send(root, nullptr, min_bounds, max_bounds, min_depth, 0, cells);
}

void Tree::cells_to_send(Cell * cell, Cell * parent, const double * min_bounds, const double * max_bounds, int min_depth, int depth, vector<Cell*> & cells){

    /* if cell is considered, we sent it */
    if(depth > min_depth){
        if(depth != min_depth + 1){
            cell->parent_idx = parent->list_idx;
        }
        else{
            cell->parent_idx = -1;
        }
        cell->list_idx = cells.size();
        cells.push_back(cell);
    }
    
    /* if subtree of cell is needed by domain, we open it */
    if(opening_criterion(cell, min_bounds, max_bounds)){
        for(int i = 0; i < 8; i++){
            if(cell->subcells[i] != nullptr){
                /* send subcells */
                cells_to_send(cell->subcells[i], cell,  min_bounds, max_bounds, min_depth, depth + 1, cells);
            }
            else{
                break;
            }
        }
    }
}
    
void Tree::delete_descendants(Cell * cell){
    /* delete subcells */
    for(int i = 0; i < 8; i++){
        if(cell->subcells[i] != nullptr){
            delete_subtree(cell->subcells[i]);
            cell->subcells[i] = nullptr;
        }
    }
}

void Tree::prune_tree(const double * min_bounds, const double * max_bounds){
    prune_tree(root, min_bounds, max_bounds);
}
    
void Tree::prune_tree(Cell * cell, const double * min_bounds, const double * max_bounds){
    /* if cell is not leaf */
    if(cell->subcells[0] != nullptr){
        /* if subtree not needed by domain delete it*/
        if(!opening_criterion(cell, min_bounds, max_bounds)){
            delete_descendants(cell);           
        }
        /* else call prune recursively */
        else{
            for(int i = 0; i < 8; i++){
                if(cell->subcells[i] != nullptr){
                    prune_tree(cell->subcells[i], min_bounds, max_bounds);
                }
            }
        }
    }
}

array<double, 3> Tree::compute_force(int bodyIndex) {
    return compute_force(root, bodyIndex);
}

array<double, 3> Tree::compute_force(const Cell* cell, int bodyIndex) {

    /* if cell is non empty and we dont want to open it or it is leaf cell */
    auto bodies = bodyManager->localBodies;
    auto pos = bodies.position.data() + bodyIndex * 4 * sizeof(double);
    if (cell->m != 0 and (!opening_criterion(cell, pos, pos) or cell->subcells[0] == nullptr)) {
        if (cell->body != bodyIndex) {
            /* evaluate force */
            return eval_force(cell->rm, cell->m, pos, bodies.mass[bodyIndex]);
        }
        /* is this really needed? */
        else {
            return { {0, 0, 0} };
        }
    }

    /* accumulate force evaluation from subcells */
    array<double, 3> ftot = { {0, 0, 0} };
    for (int i = 0; i < 8; i++) {
        if (cell->subcells[i] != nullptr) {
            /* call force computation recursively */
            array<double, 3> f = compute_force(cell->subcells[i], bodyIndex);

            /* accumulate */
            for (int c = 0; c < 3; c++) {
                ftot[c] += f[c];
            }
        }
        else {
            break;
        }
    }
    return ftot;
}
    
double Tree::mass(){
    return root->m;
}

int Tree::size(bool complete_tree){
    return size(root, complete_tree); 
}

int Tree::size(Cell * c, bool complete_tree){
    int s = 0;
    int i;
    for(i = 0; i < 8; i++){
        if(c->subcells[i] != nullptr){
            s += size(c->subcells[i], complete_tree);
        }
        else{
            break;
        }
    }
    if(i == 0 and c->m != 0){
        return 1;
    }
    else{
        if(complete_tree){
            s++;
        }
        return s; 
    }
}

bool coord_in_cell(Cell * cell, const double * pos){
    for(int c = 0; c < 3; c++){
        /* if coordinate does not lie in cell boundary */
        if(pos[c] < cell->min_bounds[c] or pos[c] > cell->max_bounds[c]){
            return false;
        }
    }
    return true;
}

bool same_cell(Cell * cell, const double * min_bounds, const double * max_bounds){
    for(int c = 0; c < 3; c++){
        /* if coordinate boundaries or not same */
        if(min_bounds[c] != cell->min_bounds[c] or max_bounds[c] != cell->max_bounds[c]){
            return false;
        }
    }
    return true;
}

double cell_volume(const Cell * cell){
    double v = 0;
    double tmp;
    for(int c = 0; c < 3; c++){
       /* coordinate side length */
       tmp = (cell->max_bounds[c] - cell->min_bounds[c]);
       if(v == 0){
            v = tmp;
       }
       else{
            v *= tmp;
       }
    }
    return v;
}


bool bounds_in_cell(Cell * cell, const double * min_bounds, const double * max_bounds){
    for(int c = 0; c < 3; c++){
        /* if coordinate boundary does not fit in cell*/
        if(cell->min_bounds[c] > min_bounds[c] or cell->max_bounds[c] < max_bounds[c]){
            return false;
        }
    }
    return true;
}
